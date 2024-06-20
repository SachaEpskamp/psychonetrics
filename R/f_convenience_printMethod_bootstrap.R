# psychonetrics print method
 
setMethod(f = "show",
 signature = "psychonetrics_bootstrap",
definition = function(object){
  
    # version:
  version <- read.dcf(file=system.file("DESCRIPTION", package="psychonetrics"),
                                  fields="Version")
  
  # Obtain some information:
  n_boots <- object@n_success
  n_fail <- object@n_fail
  boot_sub <-  object@boot_sub
  boot_resample <-  object@boot_resample
  last_run <- format(as.POSIXct(max(sapply(object@models,function(x)x@log[[length(x@log)]]@time))), "%Y-%m-%d %H:%M:%S")
  
  # Make super cool header:
  # cat(
  #   paste0("\t\t#########################\n",
  #          "\t\t## psychonetrics model ##\n",
  #          "\t\t#########################\n\n"))
  psychonetrics_print_logo()
  
  # output some general stats:
  cat("General:",
    "\n\t- Aggregated bootstrap results!",
    "\n\t- psychonetrics version:",version,
    "\n\t- Last bootstrap sample edited at:",last_run)
  
  # output some sample specific stats:
  cat("\n\nBootstrap design:",
      "\n\t- Number of included bootstraps:",n_boots,
      "\n\t- Number of removed (non-converged) bootstraps:",n_fail,
      "\n\t- Proportion of cases sampled",boot_sub,
      "\n\t- Sampling with replacement:",boot_resample)
  
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
      "\n\t- Submodel used:",submod)
  
  
    # Next steps:
    cat("\n\nTips:",
        "\n\t- Use 'psychonetrics::fit' to inspect model fit",
        "\n\t- Use 'psychonetrics::parameters' to inspect model parameters",
        "\n\t- Use 'psychonetrics::CIplot' to plot bootstrapped CIs"
    )
  # Newline to end:
  cat("\n")
  
})