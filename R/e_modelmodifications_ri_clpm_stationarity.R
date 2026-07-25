
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
    
  } else if (stationary %in% c("contemporaneous", "innovation")) {

    # Free latent-(co)variance parameter numbers ("_zeta" matrices):
    if (stationary == "contemporaneous"){
      parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row != x@parameters$col))
    } else {
      parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row == x@parameters$col))
    }

    # The latents are ordered wave-1, ..., wave-nTime (nVar latents each),
    # followed by the nVar random intercepts. Derive each parameter's
    # wave-block and within-block coordinates structurally, instead of the
    # earlier positional reshape (matrix(parNums, ncol = nTime + 1)), which
    # silently recycled and corrupted the constraints when parameters had
    # been pruned or fixed. Wave 1 is deliberately left free (its block is
    # an exogenous T1 covariance, not an innovation covariance), and the
    # random-intercept block is left free as well; only the wave-2+
    # innovation blocks are equated (identical to the old behavior for
    # unpruned models):
    parRow <- x@parameters$row[parNums]
    parCol <- x@parameters$col[parNums]
    parMat <- x@parameters$matrix[parNums]
    block_row <- ifelse(parRow > nVar * nTime, nTime + 1, ceiling(parRow / nVar))
    block_col <- ifelse(parCol > nVar * nTime, nTime + 1, ceiling(parCol / nVar))

    # Keep within-block parameters of the wave-2+ innovation blocks:
    keep <- block_row == block_col & block_row >= 2 & block_row <= nTime
    parNums <- parNums[keep]
    within_row <- parRow[keep] - (block_row[keep] - 1) * nVar
    within_col <- parCol[keep] - (block_row[keep] - 1) * nVar
    elementid <- paste0(parMat[keep], "[", within_row, ",", within_col, "]")

    # Group the same within-block element across waves:
    groups <- split(seq_along(parNums), elementid)

    # Order groups as the parameters appear in the parameter table:
    groups <- groups[order(vapply(groups, function(i) min(parNums[i]), numeric(1)))]

    # Constrain equal, requiring each element to be free in all waves:
    for (grp in seq_along(groups)){
      ind <- groups[[grp]]
      if (length(ind) != nTime - 1){
        stop("Cannot impose stationary ",
             ifelse(stationary == "innovation", "innovation variances", "contemporaneous covariances"),
             ": element ", names(groups)[grp], " (within-block indices) is free in only ",
             length(ind), " of the ", nTime - 1, " innovation blocks (waves 2-", nTime,
             "). This can happen after prune() or fixpar(); free or fix this parameter in all waves first.")
      }
      x <- parequal(x, parNums[ind])
    }

  } else stop(paste0("stationary = '",stationary,"' not implemented!"))


  # Return model:
  return(x)
}


# --------------------------------------------------------------------------
# Guarded decision-tree stationarity search for RI-CLPM / RI-CLPN models.
# --------------------------------------------------------------------------
# Starting from an (unconstrained) ri_clpm() / ri_clpn() model, the search
# tests the four stationarity constraints in the order
#
#   temporal -> contemporaneous -> intercepts -> innovation
#
# and closes with the wave-1 endogenous panel(G)VAR model. Up to psychonetrics
# 0.17.1 this was a strict chain that stopped at the first rejection. Since
# 0.17.2 the search is a *guarded decision tree* that follows the targets of
# inference:
#
#   * "temporal" is the only hard GATE. The temporal (cross-lagged) effects
#     are the primary target of inference: if they are not stationary there is
#     no single lag-1 network to report and the target of inference changes to
#     a set of per-transition temporal networks. The search then stops and the
#     unconstrained base model is selected.
#
#   * "contemporaneous", "intercepts" and "innovation" are nuisance- or
#     secondary-structure questions. A rejection is RELEASED: the constraint
#     is not imposed (that parameter set stays wave-specific) but the search
#     continues, so every one of these constraints is tested exactly once,
#     conditional only on the constraints that survived their own test. This
#     reaches endpoints the old chain could never reach (e.g. wave-varying
#     contemporaneous relations with stationary intercepts and innovation
#     variances).
#
#   * "wave1_endogenous" (the panel(G)VAR model) is gated on COVARIANCE
#     stationarity: it requires the temporal, contemporaneous and innovation
#     constraints to be retained. A wave-varying innovation covariance
#     structure has no stationary distribution, so there is nothing for an
#     endogenous first wave to be drawn from and the step is reported as "not
#     applicable". When the intercepts were released, the panel model is fit
#     with stationary_means = FALSE, so that both models carry a saturated
#     (wave-specific) mean structure and only the wave-1 exogenous block
#     differs. In both cases the panel model is the CONSTRAINED special case:
#     it has m(m+1)/2 FEWER free parameters (higher df) than the working
#     RI-CLPM, namely the freed wave-1 exogenous within-person covariance
#     block.
#
# The correlation parameterization used by ri_clpm() (type = "cor" by default)
# makes the contemporaneous (correlation) and innovation-variance (standard
# deviation) steps two cleanly separated, nested constraints; the innovation
# (co)variances are never equated with the exogenous wave-1 block nor with the
# between-person (random intercept) block.
#
# on_rejection = "stop" reproduces the pre-0.17.2 stop-at-first-rejection
# chain exactly. It is the same tree with an early break, so there is a single
# code path.
ri_clpm_search <- function(
    x,                                    # An ri_clpm() model (run or unrun)
    criterion = c("BIC", "AIC", "Chisq"), # Selection criterion
    alpha = 0.05,                         # Significance level (criterion = "Chisq")
    include_panelvar = TRUE,              # Also test the panel(G)VAR (wave-1 endogenous) model
    verbose = TRUE,
    on_rejection = c("release", "stop"),  # "release" = guarded tree, "stop" = old chain
    ...                                   # Passed to runmodel()
){
  stopifnot(is(x, "psychonetrics"))
  if (!identical(x@submodel, "RI_CLPM")){
    stop("'x' must be an ri_clpm() model.")
  }
  criterion    <- match.arg(criterion)
  on_rejection <- match.arg(on_rejection)

  # The three release-and-continue steps follow the temporal gate:
  steps    <- c("temporal", "contemporaneous", "intercepts", "innovation")
  allsteps <- c(steps, "wave1_endogenous")

  # --- Verdict bookkeeping --------------------------------------------------
  # One row per step; every step that is never reached keeps "not tested".
  verdicts <- data.frame(
    step       = allsteps,
    tested     = FALSE,
    decision   = "not tested",
    Chisq_diff = NA_real_,
    DF_diff    = NA_real_,
    p_value    = NA_real_,
    dAIC       = NA_real_,
    dBIC       = NA_real_,
    note       = "",
    stringsAsFactors = FALSE
  )
  setverdict <- function(step, ...){
    vals <- list(...)
    i <- match(step, verdicts$step)
    for (nm in names(vals)) verdicts[[nm]][i] <<- vals[[nm]]
    invisible(NULL)
  }

  # Criterion value of a model:
  critval <- function(m){
    switch(criterion,
           AIC   = m@fitmeasures$aic.ll,
           BIC   = m@fitmeasures$bic,
           Chisq = m@fitmeasures$chisq)
  }

  # Fit validity. A candidate that could not be fitted, did not converge, or
  # produced a negative / non-finite chi-square or a non-finite criterion is
  # never retained: its constraint counts as rejected. Returns "" when valid,
  # otherwise a short reason for the verdict table.
  invalidity <- function(m){
    if (is.null(m)) return("model could not be fitted")
    if (!isTRUE(m@computed)) return("model did not converge")
    fm <- m@fitmeasures
    if (!isTRUE(is.finite(fm$chisq))) return("non-finite chi-square")
    if (fm$chisq < -1e-4) return("negative chi-square")
    if (!isTRUE(is.finite(fm$df))) return("non-finite degrees of freedom")
    if (!isTRUE(is.finite(critval(m)))) return(paste0("non-finite ", criterion))
    ""
  }

  # Compare a more constrained candidate to the current working model:
  evaluate <- function(cur, cand){
    cd <- cand@fitmeasures$chisq  - cur@fitmeasures$chisq
    dd <- cand@fitmeasures$df     - cur@fitmeasures$df
    da <- cand@fitmeasures$aic.ll - cur@fitmeasures$aic.ll
    db <- cand@fitmeasures$bic    - cur@fitmeasures$bic
    p  <- if (isTRUE(dd > 0) && isTRUE(is.finite(cd))){
      pchisq(max(cd, 0), dd, lower.tail = FALSE)
    } else NA_real_
    keep <- if (criterion == "Chisq"){
      is.na(p) || p >= alpha
    } else {
      isTRUE(critval(cand) <= critval(cur))
    }
    list(keep = isTRUE(keep), Chisq_diff = cd, DF_diff = dd,
         p_value = p, dAIC = da, dBIC = db)
  }

  # An empty evaluation (used when the candidate could not be evaluated):
  noeval <- list(keep = FALSE, Chisq_diff = NA_real_, DF_diff = NA_real_,
                 p_value = NA_real_, dAIC = NA_real_, dBIC = NA_real_)

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

  # --- Base (unconstrained) model -------------------------------------------
  if (!x@computed){
    if (verbose) message("Fitting base (unconstrained) RI-CLPM ...")
    x <- runmodel(x, verbose = FALSE, ...)
  }
  models   <- list(base = x)
  current  <- x
  retained <- character(0)
  path     <- addrow(NULL, "base", "(none)", x, NA_real_, "start")

  gate_rejected <- FALSE
  stopped_early <- FALSE

  # --- The four stationarity constraints ------------------------------------
  for (st in steps){
    if (verbose) message("Testing stationary '", st, "' ...")
    cand <- tryCatch(
      runmodel(ri_clpm_stationary(current, stationary = st), verbose = FALSE, ...),
      error = function(e){
        warning("Stationarity step '", st, "' failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })

    bad <- invalidity(cand)
    if (nzchar(bad)){
      ev   <- noeval
      note <- bad
    } else {
      ev   <- evaluate(current, cand)
      note <- if (isTRUE(ev$Chisq_diff < -1e-4)){
        "negative chi-square difference (check convergence)"
      } else ""
    }

    # Record the fitted candidate (also when rejected):
    if (!is.null(cand) && isTRUE(cand@computed)){
      models[[st]] <- cand
      path <- addrow(path, st, st, cand, ev$p_value,
                     if (ev$keep) "accept" else "reject")
    }

    decision <- if (ev$keep) "retained" else if (st == "temporal") "rejected (gate)" else "released"
    setverdict(st, tested = TRUE, decision = decision,
               Chisq_diff = ev$Chisq_diff, DF_diff = ev$DF_diff,
               p_value = ev$p_value, dAIC = ev$dAIC, dBIC = ev$dBIC,
               note = note)

    if (ev$keep){
      current  <- cand
      retained <- c(retained, st)
    } else if (st == "temporal"){
      # The gate: the primary target of inference changes, stop the search.
      if (verbose) message("Stationary temporal effects rejected: this is a gate, stopping the search.")
      gate_rejected <- TRUE
      break
    } else if (on_rejection == "stop"){
      stopped_early <- TRUE
      break
    }
  }

  # --- Step "wave1_endogenous": the panel(G)VAR model -----------------------
  # Gated on covariance stationarity of the within-person process.
  nVar            <- nrow(x@extramatrices$vars)
  expected_df_gap <- nVar * (nVar + 1) / 2
  free_means      <- !("intercepts" %in% retained)
  panel_selected  <- FALSE

  cov_stationary <- all(c("temporal", "contemporaneous", "innovation") %in% retained)

  if (gate_rejected || stopped_early){
    # Leave as "not tested" (the search already stopped).
  } else if (!include_panelvar){
    setverdict("wave1_endogenous", tested = FALSE, decision = "not tested",
               note = "include_panelvar = FALSE")
  } else if (!cov_stationary){
    setverdict("wave1_endogenous", tested = FALSE, decision = "not applicable",
               note = paste0("the within-person covariance process is not stationary",
                             " (contemporaneous and/or innovation structure released),",
                             " so the endogenous wave-1 model is undefined"))
  } else {
    if (verbose){
      message("Testing wave-1 endogenous panel(G)VAR",
              if (free_means) " (wave-specific means)" else "", " ...")
    }
    type <- x@extramatrices$ri_clpm_type
    if (is.null(type)) type <- x@types$latent
    design <- x@extramatrices$vars

    pv <- tryCatch(
      withCallingHandlers(
        runmodel(panelvar(vars = design, sampleStats = x@sample,
                          within_latent = type, between_latent = type,
                          baseline = "none",
                          stationary_means = !free_means),
                 verbose = FALSE),
        warning = function(w){
          # The panel(G)VAR step uses no baseline (incremental fit not needed here):
          if (grepl("No baseline model found", conditionMessage(w))) invokeRestart("muffleWarning")
        }),
      error = function(e){
        warning("panelVAR step failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })

    bad <- invalidity(pv)
    if (nzchar(bad)){
      ev   <- noeval
      note <- bad
    } else {
      ev   <- evaluate(current, pv)
      note <- ""

      # The panel model is the constrained special case: it should have
      # exactly m(m+1)/2 MORE degrees of freedom (the freed wave-1 exogenous
      # within-person covariance block). Anything else means the constraint
      # bookkeeping of the working model is not what this step assumes.
      if (!isTRUE(all.equal(ev$DF_diff, expected_df_gap))){
        warning("The wave-1 (panel(G)VAR) comparison has a degrees-of-freedom gap of ",
                ev$DF_diff, " instead of the expected m(m+1)/2 = ", expected_df_gap,
                ". The panel model should be the constrained special case of the working ",
                "RI-CLPM, differing only in the wave-1 exogenous covariance block. ",
                "This indicates that the constraint bookkeeping went wrong upstream ",
                "(for example after prune() or fixpar()); interpret this step with care.",
                call. = FALSE)
        note <- paste0("unexpected df gap (", ev$DF_diff, " instead of ", expected_df_gap, ")")
      }
      if (isTRUE(ev$Chisq_diff < -1e-4)){
        note <- if (nzchar(note)){
          paste0(note, "; negative chi-square difference (check convergence)")
        } else "negative chi-square difference (check convergence)"
      }
    }

    if (!is.null(pv) && isTRUE(pv@computed)){
      models[["panelVAR"]] <- pv
      path <- addrow(path, "panelVAR", "wave1_endogenous", pv, ev$p_value,
                     if (ev$keep) "accept" else "reject")
    }

    setverdict("wave1_endogenous", tested = TRUE,
               decision = if (ev$keep) "retained" else "released",
               Chisq_diff = ev$Chisq_diff, DF_diff = ev$DF_diff,
               p_value = ev$p_value, dAIC = ev$dAIC, dBIC = ev$dBIC,
               note = note)

    if (ev$keep){
      current        <- pv
      panel_selected <- TRUE
    }
  }

  # --- Name of the selected model -------------------------------------------
  selected_name <- if (panel_selected){
    paste0("panelVAR", if (free_means) " (means free)" else "")
  } else if (length(retained) == 0){
    "base"
  } else {
    paste0(paste(retained, collapse = "+"),
           if (free_means) " (means free)" else "")
  }

  # Full comparison table. The models are sorted by df; when constraints were
  # released, consecutive rows need not be nested, so the (informative)
  # negative chi-square difference warning of compare() is muffled here:
  comparison <- tryCatch(
    withCallingHandlers(
      do.call(compare, models),
      warning = function(w){
        if (grepl("Negative chi-square difference", conditionMessage(w))) invokeRestart("muffleWarning")
      }),
    error = function(e) NULL)

  out <- list(
    selected      = current,
    selected_name = selected_name,
    models        = models,
    comparison    = comparison,
    path          = path,
    verdicts      = verdicts,
    criterion     = criterion,
    alpha         = alpha,
    on_rejection  = on_rejection
  )
  class(out) <- "ri_clpm_search"
  if (verbose) message("Selected model: ", selected_name)
  out
}


# Print method for ri_clpm_search objects:
print.ri_clpm_search <- function(x, ...){
  # Small formatting helpers (blank instead of NA):
  num <- function(v, d = 2) ifelse(is.na(v), "", formatC(v, format = "f", digits = d))

  cat("RI-CLPM / RI-CLPN stationarity search\n")
  cat("Criterion: ", x$criterion,
      if (x$criterion == "Chisq") paste0(" (alpha = ", x$alpha, ")") else "",
      "\n", sep = "")
  if (!is.null(x$on_rejection)){
    cat("On rejection: ", x$on_rejection,
        if (identical(x$on_rejection, "release")){
          " (guarded tree: only the temporal gate stops the search)"
        } else " (stop at first rejection)", "\n", sep = "")
  }

  cat("\nModels fitted:\n")
  pp <- x$path
  pp$AIC    <- round(pp$AIC, 2)
  pp$BIC    <- round(pp$BIC, 2)
  pp$Chisq  <- round(pp$Chisq, 2)
  pp$p_diff <- ifelse(is.na(pp$p_diff), "", formatC(pp$p_diff, format = "f", digits = 3))
  print(pp, row.names = FALSE)

  if (!is.null(x$verdicts)){
    cat("\nVerdicts:\n")
    vv <- x$verdicts
    notes <- vv$note
    vv$Chisq_diff <- num(vv$Chisq_diff)
    vv$DF_diff    <- ifelse(is.na(vv$DF_diff), "", as.character(vv$DF_diff))
    vv$p_value    <- ifelse(is.na(vv$p_value), "", formatC(vv$p_value, format = "f", digits = 3))
    vv$dAIC       <- num(vv$dAIC)
    vv$dBIC       <- num(vv$dBIC)
    vv$note       <- NULL
    print(vv, row.names = FALSE)

    if (any(nzchar(notes))){
      cat("\nNotes:\n")
      for (i in which(nzchar(notes))){
        cat(" - ", x$verdicts$step[i], ": ", notes[i], "\n", sep = "")
      }
    }

    # The temporal gate: explain why the search stopped.
    gate <- x$verdicts$decision[x$verdicts$step == "temporal"]
    if (length(gate) == 1 && identical(gate, "rejected (gate)")){
      cat("\nThe temporal (cross-lagged) effects are NOT stationary.\n")
      cat("This step is a gate, so the search stopped and the unconstrained model is selected.\n")
      cat("Why: the temporal effects are the primary target of inference. If they differ over\n")
      cat("waves there is no single lag-1 temporal network to report, and the target of\n")
      cat("inference changes to a set of per-transition temporal networks (one per wave-to-\n")
      cat("wave step). The remaining constraints were not tested, because each of them would\n")
      cat("only be interpretable conditional on a temporal structure that does not hold.\n")
      cat("Before concluding non-stationarity, check whether the waves are EQUALLY SPACED in\n")
      cat("time: unequal measurement intervals produce wave-varying lag-1 effects even when\n")
      cat("the underlying process is perfectly stationary (a longer interval implies weaker\n")
      cat("lag-1 effects), and so mimic temporal non-stationarity.\n")
    }
  }

  cat("\nSelected model: ", x$selected_name, "\n", sep = "")
  invisible(x)
}
