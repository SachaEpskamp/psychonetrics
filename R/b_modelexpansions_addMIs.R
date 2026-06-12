# Full function:
addMIs <- function(x, matrices = "all", type =  c("normal","free","equal"),verbose,analyticFisher=TRUE){
  if (missing(verbose)){
    verbose <- x@verbose
  }

  # Check the matrices argument:
  if (!identical(matrices, "all")){
    if (!all(matrices %in% x@parameters$matrix)){
      stop(paste0("The following matrices are not part of the model: ",
                  paste0("'", setdiff(matrices, x@parameters$matrix), "'", collapse = ", ")))
    }
  }

  # full <- TRUE
  # if (full){
  if (verbose){
    message("Computing modification indices...")
  }

  tryres <- try({

    if ("normal" %in% type){

      x <-  x %>% addMIs_inner_full(matrices = matrices, type = "normal",analyticFisher=analyticFisher)
    }

    if (nrow(x@sample@groups) > 1){
      if ("free" %in% type){
        # if (verbose){
        #   message("Computing constrain-free modification indices...")
        # }
        x <- x %>% addMIs_inner_full(matrices = matrices, type = "free")
      }
      if ("equal" %in% type){

        # if (verbose){
        #   message("Computing group-constrained indices...")
        # }
        x <- x %>% addMIs_inner_full(matrices = matrices, type = "equal")
      }

    }

  })
  
  if (is(tryres, "try-error")){
    stop("Failed to compute modification indices. Try a different optimizer with setoptimizer(...) or report the problem on github.com/SachaEpskamp/psychonetrics.")
  }
  
  return(x)
  
}


# Add the modification indices (FULL VERSION):
addMIs_inner_full <- function(x, matrices = "all", type =  c("normal","free","equal"),analyticFisher=TRUE){
  type <- match.arg(type)

  # Which matrices to include in the MIs:
  useMatrices <- if (identical(matrices, "all")) unique(x@parameters$matrix) else matrices
  
  # If no constrained parameters, nothing to do!
  if (!any(x@parameters$par == 0) & !any(duplicated(x@parameters$par))){
    return(x)
  }
  
  # Clear old MIs:
  if (type == "normal"){
    x@parameters$mi[] <- 0
    x@parameters$pmi[] <- NA
    x@parameters$epc[] <- NA
  } else if (type == "free"){
    x@parameters$mi_free[] <- 0
    x@parameters$pmi_free[] <- NA
    x@parameters$epc_free[] <- NA
    # Joint score-test columns (added in 0.15.4). Create if absent for older models:
    if (is.null(x@parameters$mi_free_joint)) x@parameters$mi_free_joint <- NA_real_
    if (is.null(x@parameters$pmi_free_joint)) x@parameters$pmi_free_joint <- NA_real_
    if (is.null(x@parameters$df_free_joint)) x@parameters$df_free_joint <- NA_real_
    x@parameters$mi_free_joint[] <- NA_real_
    x@parameters$pmi_free_joint[] <- NA_real_
    x@parameters$df_free_joint[] <- NA_real_
  } else {
    x@parameters$mi_equal[] <- 0
    x@parameters$pmi_equal[] <- NA
    x@parameters$epc_equal[] <- NA
  }
  
  # Sample size:
  n <- sum(x@sample@groups$nobs)
  
  # Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but constrained per group #
  # Fully free:
  # Copy the model:
  modCopy <- x
  
  # Obtain the full set of parameters that are constrained across all groups:
  if (type == "equal"){
    sum <- modCopy@parameters %>% filter(.data[['matrix']] %in% useMatrices) %>%
      group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% summarize(anyConstrained = any(.data[['fixed']])) %>%
      filter(drop(.data[['anyConstrained']]))
    # Add a unique number to each:
    sum$par2 <-  max(modCopy@parameters$par) + seq_len(nrow(sum))

    # Left join back (rows of excluded matrices have NA par2 and keep their par):
    modCopy@parameters <- modCopy@parameters %>% left_join(sum,by=c("matrix","row","col")) %>%
      mutate(par = ifelse(.data[['identified']],0,ifelse(.data[['par']]==0 & !is.na(.data[['par2']]),.data[['par2']],.data[['par']])))

  } else {
    # Add free parameter numbers (only in the requested matrices):
    toFree <- modCopy@parameters$par==0 & !modCopy@parameters$identified & modCopy@parameters$matrix %in% useMatrices
    modCopy@parameters$par[toFree] <- max(modCopy@parameters$par) + seq_len(sum(toFree))

    # For each group, free all parameters from equality constraints:
    if (type == "free"){
      if (nrow(modCopy@sample@groups)>1){
        for (g in 2:nrow(modCopy@sample@groups)){
          toFree <- modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified & modCopy@parameters$matrix %in% useMatrices
          modCopy@parameters$par[toFree] <- max(modCopy@parameters$par) + seq_len(sum(toFree))
        }
      }
    }

  }

  # If the matrices filter left nothing to free, there is nothing to test:
  if (max(modCopy@parameters$par) == max(x@parameters$par)){
    return(x)
  }
  
  # Identify:
  # modCopy <- identify(modCopy)
  
  # Remake the model matrix:
  # modCopy@extramatrices$M <- Mmatrix(modCopy@parameters)
  # 
  if (modCopy@cpp){
    modCopy@extramatrices$M   <- Mmatrix_cpp(modCopy@parameters )
  } else {
    modCopy@extramatrices$M  <- Mmatrix(modCopy@parameters)  
  }
  
  
  # Check if gradient and hessian are present:
  # gradient <- !is.null(x@fitfunctions$gradient)
  # hessian <- !is.null(x@fitfunctions$hessian)
  # information <- !is.null(x@fitfunctions$information)
  
  # Compute a gradient:
  if (modCopy@cpp){
    g <- psychonetrics_gradient_cpp(parVector(modCopy), modCopy)
  } else {
    g <- psychonetrics_gradient(parVector(modCopy), modCopy)
    
  }
  # if (gradient){
  #   g <- x@fitfunctions$gradient(parVector(modCopy), modCopy)
  # } else {
  #   g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  
  # Compute the Fisher information:
  # if (information){
  #   H <- x@fitfunctions$information(modCopy)
  # } else if (hessian){
  #   H <- x@fitfunctions$hessian(parVector(modCopy), modCopy)
  # } else if (gradient){
  #   H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(modCopy), model=modCopy) 
  # } else {
  #   H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  # Total sample size:
  nTotal <- sum(x@sample@groups$nobs)
  
  # FIXME: 4 * n could be nicer here probably
  if (!analyticFisher){
    # The C++ analytic path does not honour analyticFisher (its 2nd argument is
    # useM, not analytic), so route to the numeric R implementation here:
    H <- 4 * nTotal * as(numeric_FisherInformation(modCopy), "matrix")
  } else if (modCopy@cpp){
    H <- 4 * nTotal * as(psychonetrics_FisherInformation_cpp(modCopy, analyticFisher), "matrix")
  } else {
    H <- 4 * nTotal * as(psychonetrics_FisherInformation(modCopy, analyticFisher), "matrix")
  }
  
  
  # For every new parameter:
  curMax <- max(x@parameters$par)
  
  ### NEW
  curInds <- seq_len(curMax)
  newInds <- curMax + seq_len(max(modCopy@parameters$par) - curMax)
  V <-  H[newInds,newInds] - H[newInds,curInds,drop=FALSE] %*% solve_symmetric(H[curInds,curInds]) %*% H[curInds,newInds,drop=FALSE]
  V.diag <- diag(V)
  idx <- which(V.diag < sqrt(.Machine$double.eps))
  if(length(idx) > 0L) {
    V.diag[idx] <- as.numeric(NA)
  }
  # How many in total?
  nTotalPars <- length(c(curInds,newInds))
  
  # 
  # 
  # # Effective N:
  # if (nrow(x@sample@groups) > 1){
  #   par <- modCopy@parameters
  #   # par$id[!modCopy@parameters$identified] <- seq_len(sum(!modCopy@parameters$identified))
  #   par <- par %>% left_join(modCopy@sample@groups, by = c("group_id" = "id")) %>%
  #     group_by(par) %>% summarize(Neff = sum(nobs))
  #   Neff <- numeric(max(par$par))
  #   Neff[par$par[par$par!=0]] <- par$Neff[par$par!=0]
  # } else {
  #   Neff <- rep(x@sample@groups$nobs[1],nTotalPars)
  # }
  
  
  # MIs:
  # All MIs:
  mi <- numeric(nTotalPars)
  
  # mi[newInds] <- ifelse(abs(V.diag) < 1e-10,0,((-Neff[newInds]*g[newInds])^2)/V.diag)
  mi[newInds] <- ifelse(abs(V.diag) < 1e-10,0,((-nTotal*g[newInds])^2)/V.diag)
  if (length(curInds) > 0){
    # mi[curInds] <- ((-Neff[curInds]*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
    mi[curInds] <- ((-nTotal*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
  }
  p <- pchisq(mi,df = 1,lower.tail = FALSE)     
  
  # Compute epc:
  # d <- 0.5 * (-1 * Neff) * g
  d <- 0.5 * (-1 * nTotal) * g
  # needed? probably not; just in case
  d[which(abs(d) < 1e-15)] <- 1.0
  
  # Expected parameter change:
  epc <-   mi / d 
  
  # Which to fill (restricted to the requested matrices; no-op for "all"):
  fillInds <- match(c(curInds,newInds),modCopy@parameters$par)
  fillSel <- !is.na(fillInds) & modCopy@parameters$matrix[fillInds] %in% useMatrices
  if (type == "normal"){
    x@parameters$mi[fillInds[fillSel]] <- round(mi[fillSel],10) # round(mi, 3)
    x@parameters$pmi[fillInds[fillSel]] <- round(p[fillSel],10)
    x@parameters$epc[fillInds[fillSel]] <- round(epc[fillSel],10)
  } else if (type == "free"){
    x@parameters$mi_free[fillInds[fillSel]] <- round(mi[fillSel],10) # round(mi, 3)
    x@parameters$pmi_free[fillInds[fillSel]] <- round(p[fillSel],10)
    x@parameters$epc_free[fillInds[fillSel]] <- round(epc[fillSel],10)

    # --- Joint score test (Lagrange Multiplier) for releasing each equality constraint ---
    # We delegate to .equalityScoreTestInner (in f_convenience_equalityScoreTest.R) which
    # does its own clean augmentation (only releasing equality constraints, not freeing
    # zero pars). The result populates mi_free_joint / pmi_free_joint / df_free_joint
    # columns; the original per-parameter mi_free values above are left as-is.
    # joint_only = TRUE: addMIs only writes the joint columns (mi_free_joint
    # etc.); the univariate score-test rows are not consumed here. Skipping
    # them avoids K * (G - 1) MASS::ginv calls per runmodel(), which dominate
    # the cost of partialprune() on multi-group models with many equality
    # constraints. The user-facing equalityScoreTest() and MIs() compute both
    # tables on demand.
    est_res <- tryCatch(.equalityScoreTestInner(x, analyticFisher = analyticFisher,
                                                method = "jacobian",
                                                joint_only = TRUE),
                        error = function(e) NULL)
    if (!is.null(est_res) && !is.null(est_res$total) && nrow(est_res$total) > 0){
      tot <- est_res$total
      # Respect the matrices filter:
      tot <- tot[tot$matrix %in% useMatrices, , drop = FALSE]
      for (i in seq_len(nrow(tot))){
        write_rows <- which(x@parameters$matrix == tot$matrix[i] &
                              x@parameters$row == tot$row[i] &
                              x@parameters$col == tot$col[i])
        x@parameters$mi_free_joint[write_rows] <- round(tot$X2[i], 10)
        x@parameters$pmi_free_joint[write_rows] <- round(tot$p.value[i], 10)
        x@parameters$df_free_joint[write_rows] <- tot$df[i]
      }
    }
  } else {
    x@parameters$mi_equal[fillInds[fillSel]] <- round(mi[fillSel],10) # round(mi,3)
    x@parameters$pmi_equal[fillInds[fillSel]] <- round(p[fillSel], 10)
    x@parameters$epc_equal[fillInds[fillSel]] <- round(epc[fillSel],10)
  }
  
  return(x)
}