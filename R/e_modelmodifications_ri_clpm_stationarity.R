
# FUnction to make the RI CLPM model more stationary:
ri_clpm_stationary <- function(
    x,
    stationary = c( "intercepts", "contemporaneous", "innovation", "temporal")
){
  stopifnot(is(x, "psychonetrics"))
  
  # Check stationary 
  stationary <- match.arg(stationary)
  
  # Verbose:
  verbose <- x@verbose
  
  # Extra matrices:
  vars <- x@extramatrices$vars
  varsDF <- x@extramatrices$varsDF
  
  # Number of variables:
  nVar <- nrow(vars)
  
  # Number of measurements:
  nTime <- ncol(vars)
  
  # Number of total observed:
  n_observed <- nVar * nTime
  
  # Set stationary:
  if (stationary == "intercepts"){
    
    # Average intercept:
    avg_int <- rowMeans(matrix(getmatrix(x, "nu", group = 1),nrow=nVar))
    
    # Fix all nu parameters:
    x <- fixpar(x, matrix = "nu", row = seq_len(n_observed), col = 1)
    
    # Free nu_eta parameters for random intercepts:
    x <- freepar(x, matrix = "nu_eta", row = n_observed + seq_len(nVar), col = 1)
    
    # Set starting values:
    x@parameters$est[x@parameters$matrix == "nu_eta" & !x@parameters$fixed] <- avg_int
    
  } else if (stationary == "temporal") {
    
    # Beta parameter numbers:
    parNums <- which(x@parameters$matrix == "beta" & !x@parameters$fixed)
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar^2){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime-1)
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  } else if (stationary == "contemporaneous") {

    # Beta parameter numbers:
    parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row != x@parameters$col))
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar*(nVar-1)/2){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime+1)
    
    # cut out time = 1 and RI:
    parNums <- parNums[,-c(1,ncol(parNums)),drop=FALSE]
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  }  else if (stationary == "innovation") {
    
    # Beta parameter numbers:
    parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row == x@parameters$col))
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime+1)
    
    # cut out time = 1 and RI:
    parNums <- parNums[,-c(1,ncol(parNums)),drop=FALSE]
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  } else stop(paste0("stationary = '",stationary,"' not implemented!"))


  # Return model:
  return(x)
}


# Sequential stationarity search for RI-CLPM models.
# Starts from an (unconstrained) ri_clpm model and adds, one at a time,
# stationarity constraints on (1) intercepts, (2) temporal effects,
# (3) contemporaneous relations, (4) innovation variances, and finally
# (5) the panelVAR model (wave 1 treated as the stationary distribution).
# At each step the more constrained model is compared to the currently
# selected one and the constraint is retained only if it does not worsen
# fit according to 'criterion'. The search stops at the first rejection.
ri_clpm_search <- function(
    x,                                    # An ri_clpm() model (run or unrun)
    criterion = c("BIC", "AIC", "Chisq"), # Selection criterion
    alpha = 0.05,                         # Significance level (criterion = "Chisq")
    include_panelvar = TRUE,              # Also test the panelVAR (wave-1 endogenous) model
    verbose = TRUE,
    ...                                   # Passed to runmodel()
){
  stopifnot(is(x, "psychonetrics"))
  if (!identical(x@submodel, "RI_CLPM")){
    stop("'x' must be an ri_clpm() model.")
  }
  criterion <- match.arg(criterion)

  # Criterion value of a model:
  critval <- function(m){
    switch(criterion,
           AIC   = m@fitmeasures$aic.ll,
           BIC   = m@fitmeasures$bic,
           Chisq = m@fitmeasures$chisq)
  }

  # Decide whether to accept the more constrained candidate over current:
  decide <- function(cur, cand){
    cd <- cand@fitmeasures$chisq - cur@fitmeasures$chisq
    dd <- cand@fitmeasures$df    - cur@fitmeasures$df
    p  <- if (dd > 0) pchisq(cd, dd, lower.tail = FALSE) else NA_real_
    acc <- if (criterion == "Chisq"){
      is.na(p) || p >= alpha
    } else {
      critval(cand) <= critval(cur)
    }
    list(accept = isTRUE(acc), p = p)
  }

  # Append a row to the path table:
  addrow <- function(tab, name, constraint, m, p, decision){
    rbind(tab, data.frame(
      model      = name,
      constraint = constraint,
      DF         = m@fitmeasures$df,
      AIC        = m@fitmeasures$aic.ll,
      BIC        = m@fitmeasures$bic,
      Chisq      = m@fitmeasures$chisq,
      p_diff     = p,
      decision   = decision,
      stringsAsFactors = FALSE
    ))
  }

  # --- Base (unconstrained) model ---
  if (!x@computed){
    if (verbose) message("Fitting base (unconstrained) RI-CLPM ...")
    x <- runmodel(x, verbose = FALSE, ...)
  }
  models   <- list(base = x)
  current  <- x
  selected_name <- "base"
  path <- addrow(NULL, "base", "(none)", x, NA_real_, "start")

  # --- Sequential stationarity constraints ---
  steps <- c("intercepts", "temporal", "contemporaneous", "innovation")
  for (st in steps){
    if (verbose) message("Testing stationary '", st, "' ...")
    cand <- tryCatch(
      runmodel(ri_clpm_stationary(current, stationary = st), verbose = FALSE, ...),
      error = function(e){
        warning("Stationarity step '", st, "' failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (is.null(cand)) break
    d <- decide(current, cand)
    models[[st]] <- cand
    path <- addrow(path, st, st, cand, d$p, if (d$accept) "accept" else "reject")
    if (d$accept){
      current <- cand
      selected_name <- st
    } else {
      break
    }
  }

  # --- panelVAR step (only if fully stationary so far) ---
  if (include_panelvar && selected_name == "innovation"){
    if (verbose) message("Testing panelVAR (wave-1 endogenous) ...")
    type   <- x@extramatrices$ri_clpm_type
    if (is.null(type)) type <- x@types$latent
    design <- x@extramatrices$vars
    pv <- tryCatch(
      withCallingHandlers(
        runmodel(panelvar(vars = design, sampleStats = x@sample,
                          within_latent = type, between_latent = type,
                          baseline = "none"),
                 verbose = FALSE),
        warning = function(w){
          # The panelVAR step uses no baseline (incremental fit not needed here):
          if (grepl("No baseline model found", conditionMessage(w))) invokeRestart("muffleWarning")
        }),
      error = function(e){
        warning("panelVAR step failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(pv)){
      d <- decide(current, pv)
      models[["panelVAR"]] <- pv
      path <- addrow(path, "panelVAR", "wave1_endogenous", pv, d$p, if (d$accept) "accept" else "reject")
      if (d$accept){
        current <- pv
        selected_name <- "panelVAR"
      }
    }
  }

  # Full comparison table:
  comparison <- tryCatch(do.call(compare, models), error = function(e) NULL)

  out <- list(
    selected      = current,
    selected_name = selected_name,
    models        = models,
    comparison    = comparison,
    path          = path,
    criterion     = criterion,
    alpha         = alpha
  )
  class(out) <- "ri_clpm_search"
  if (verbose) message("Selected model: ", selected_name)
  out
}


# Print method for ri_clpm_search objects:
print.ri_clpm_search <- function(x, ...){
  cat("RI-CLPM stationarity search\n")
  cat("Criterion: ", x$criterion,
      if (x$criterion == "Chisq") paste0(" (alpha = ", x$alpha, ")") else "",
      "\n\n", sep = "")
  pp <- x$path
  pp$AIC   <- round(pp$AIC, 2)
  pp$BIC   <- round(pp$BIC, 2)
  pp$Chisq <- round(pp$Chisq, 2)
  pp$p_diff <- ifelse(is.na(pp$p_diff), "", formatC(pp$p_diff, format = "f", digits = 3))
  print(pp, row.names = FALSE)
  cat("\nSelected model: ", x$selected_name, "\n", sep = "")
  invisible(x)
}