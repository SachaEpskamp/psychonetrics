# psychonetrics print method
 
setMethod(f = "show",
 signature = "psychonetrics",
definition = function(object){
  
  # Bootstrap warning:
  if (object@sample@bootstrap){
    boot_warning()
  }
  
    # version:
  version <- read.dcf(file=system.file("DESCRIPTION", package="psychonetrics"),
                                  fields="Version")
  
  # Make super cool header:
  # cat(
  #   paste0("\t\t#########################\n",
  #          "\t\t## psychonetrics model ##\n",
  #          "\t\t#########################\n\n"))
  psychonetrics_print_logo()
  
  # output some general stats:
  cat("General:",
    "\n\t- psychonetrics version:",version,
    "\n\t- Model last edited at:",as.character(object@log[[length(object@log)]]@time))
  
  # output some sample specific stats:
  cat("\n\nSample:",
      "\n\t- Number of cases:",sum(object@sample@groups$nobs),
      "\n\t- Number of groups:",nrow(object@sample@groups),
      "\n\t- Number of observed summary statistics:",object@sample@nobs)
  
  # output some model specific stats:

  mod <- switch(
    object@model,
    # "lnm" = "Latent Network Model (LNM)",
    # "ggm" = "Gaussian graphical model (GGM)",
    # "rnm" = "Residual network model (RNM)",
    "gvar" = "Graphical vector-autoregression (GVAR)",
    "varcov" = "Variance-covariance matrix (varcov)",
    # "cholesky" = "Cholesky decomposition (cholesky)",
    "lvm" = "Latent variable model (LVM)",
    "var1" = "Lag-1 vector-autoregression (VAR1)",
    "panelvar1" = "Lag-1 panel vector auto-regression (panelvar1)",
    "dlvm1" = "Lag-1 dynamic latent variable model for panel data (dlvm1)",
    "tsdlvm1" = "Lag-1 dynamic latent variable model for time-series data (tsdlvm1)",
    "meta_varcov" = "Variance-covariance matrix meta analysis",
    "Ising" = "Ising model",
    "ml_lvm" = "Multi-level latent variable model"
  )
  
  submod <- switch(
    object@submodel,
    "none" = "none",
    "lnm" = "Latent Network Model (LNM)",
    "ggm" = "Gaussian graphical model (GGM)",
    "rnm" = "Residual network model (RNM)",
    "gvar" = "Graphical vector-autoregression (GVAR)",
    "cholesky" = "Cholesky decomposition (cholesky)",
    "sem" = "Structural equation model (SEM)",
    "lrnm" = "Latent & residual network model (LRNM)",
    "gvar" = "Graphical vector-autoregression (GVAR)",
    "var" = "Vector-autoregression (VAR)",
    "ml_lnm" = "Multi-level latent network model",
    "ml_rnm" = "Multi-level residual network model"
  )
  if (is.null(submod)){
    submod <- object@submodel
  }
  
  cat("\n\nModel:",
      "\n\t- Model used:",mod,
      "\n\t- Submodel used:",submod,
      "\n\t- Number of parameters:",max(object@parameters$par))
  
  # Not computeD:
  if (!object@computed){
    cat("\n\nModel has not yet been computed. Use 'runmodel' to compute parameters and fit measures.")
  } else {
    estimator <- switch(object@estimator,
            "ML" = "Maximum likelihood estimation (ML)",
            "FIML" = "Full information maximum likelihood (FIML)",
            "ULS" = "Unweighted least squares (ULS)",
            "WLS" = "Weighted least squares (WLS)",
            "DWLS" = "Diagonally weighted least squares (DWLS)",
            "PML" = "Penalized maximum likelihood estimation (PML)")
    
    # output some optimizer results:
    cat("\n\nEstimation:",
        "\n\t- Optimizer used:",object@optim$optimizer,
        "\n\t- Estimator used:",estimator,
        # "\n\t- Number of iterations:",object@optim$iterations,
        "\n\t- Message:",object@optim$message
        )

    # PML penalty details:
    if (object@estimator == "PML") {
      pen_type <- if (object@penalty$alpha == 1) "LASSO (L1)"
                  else if (object@penalty$alpha == 0) "Ridge (L2)"
                  else paste0("Elastic net (alpha = ", object@penalty$alpha, ")")
      n_penalized <- sum(object@parameters$penalty_lambda > 0 & !object@parameters$fixed)
      n_zero <- sum(object@parameters$penalty_lambda > 0 & !object@parameters$fixed & abs(object@parameters$est) < 1e-8)
      cat("\n\nPenalization:",
          "\n\t- Penalty type:", pen_type,
          "\n\t- Lambda:", object@penalty$lambda,
          "\n\t- Penalized parameters:", n_penalized,
          "\n\t- Parameters at zero:", n_zero)
    }

    # output some fit measures (inspired by Lavaan):
    if (object@estimator == "PML") {
      # Fit measures are not computed for PML (penalized objective is not a proper likelihood)
    } else if (!is.null(object@fitmeasures)){
      cat("\n\nFit:",
          "\n\t- Model Fit Test Statistic:",goodNum(object@fitmeasures$chisq),
          "\n\t- Degrees of freedom:",object@fitmeasures$df,
          "\n\t- p-value (Chi-square):",goodNum(object@fitmeasures$pvalue)
          # "\n\t- RMSEA:",goodNum(object@fitmeasures$rmsea)
      )
    } else {
      cat("\n\nFit has not yet been computed. Use 'addfit' to compute fit measures.")
    }
    # Next steps:
    if (object@estimator == "PML") {
      cat("\n\nTips:",
          "\n\t- Use 'psychonetrics::refit' for post-selection inference with standard errors",
          "\n\t- Use 'psychonetrics::parameters' to inspect model parameters",
          "\n\t- Use 'psychonetrics::penaltyVector' to inspect per-parameter penalties"
      )
    } else {
      cat("\n\nTips:",
          "\n\t- Use 'psychonetrics::compare' to compare psychonetrics models",
          "\n\t- Use 'psychonetrics::fit' to inspect model fit",
          "\n\t- Use 'psychonetrics::parameters' to inspect model parameters",
          "\n\t- Use 'psychonetrics::MIs' to inspect modification indices"
      )
    }
  }
  # Newline to end:
  cat("\n")
  
})