# psychonetrics print method
 
setMethod(f = "show",
 signature = "psychonetrics",
definition = function(object){
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
    "lnm" = "Latent Network Model (LNM)",
    "ggm" = "Gaussian graphical model (GGM)"
  )
  cat("\n\nModel:",
      "\n\t- model used:",mod,
      "\n\t- Number of parameters:",max(object@parameters$par))
  
  # Not computeD:
  if (!object@computed){
    cat("\n\nModel has not yet been computed. Use 'runmodel' to compute parameters and fit measures.")
  } else {
    estimator <- switch(object@estimator,
            "ML" = "Maximum likelihood estimation (ML)",
            "ULS" = "Unweighted least squares (ULS)",
            "WLS" = "Weighted least squares (WLS)",
            "DWLS" = "Diagonally weighted least squares (DWLS)")
    
    # output some optimizer results:
    cat("\n\nEstimation:",
        "\n\t- Optimizer used:",object@optim$optimizer,
        "\n\t- Estimator used:",estimator,
        # "\n\t- Number of iterations:",object@optim$iterations,
        "\n\t- Message:",object@optim$message
        )
    
    # output some fit measures (inspired by Lavaan):
    if (!is.null(object@fitmeasures)){
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
    cat("\n\nTips:",
        "\n\t- Use 'psychonetrics::compare' to compare psychonetrics models",
        "\n\t- Use 'psychonetrics::fit' to inspect model fit",
        "\n\t- Use 'psychonetrics::parameters' to inspect model parameters",
        "\n\t- Use 'psychonetrics::MIs' to inspect modification indices"
    )
  }
  
})