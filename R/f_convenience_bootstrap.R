# Bootstrap the data:
bootstrap <- function(x, replacement = TRUE, proportion = 1, verbose = TRUE, storedata = FALSE, baseline_saturated = TRUE){
  if (x@distribution != "Gaussian"){
    stop("Only Gaussian distirbution models supported.")
  }
  
  # Set computed to false:
  x@computed <- FALSE
  
  # Obtain the raw data:
  if (is.null(x@sample@rawdata) || nrow(x@sample@rawdata) == 0){
    stop("Raw data is needed. Please run input with 'storedata = TRUE'")
  }
  data <- x@sample@rawdata
  
  # bootstrap the data:
  nCase <- nrow(data)

  data <- data[sample(seq_len(nCase),round(proportion * nCase),replacement),]
  
  # New data:
  x@sample <- samplestats(
    data = data, # Dataset
    vars = attr(data, "vars"), # character indicating the variables Extracted if missing from data - group variable
    groups = attr(data, "groups"), # ignored if missing. Can be character indicating groupvar, or vector with names of groups
    # covs, # alternative covs (array nvar * nvar * ngroup)
    # means, # alternative means (matrix nvar * ngroup)
    # nobs, # Alternative if data is missing (length ngroup)
    missing =  attr(data, "missing"),
    fimldata = x@estimator == "FIML",
    verbose = verbose,
    storedata = storedata
  )
  
  # Number of observations:
  nVar <- nrow(x@sample@variables)
  nGroup <- nrow(x@sample@groups)
  x@sample@nobs <-  
    nVar * (nVar+1) / 2 * nGroup + # Covariances per group
    nVar * nGroup
  
  # Reset start values:
  x <- emergencystart(x)
  
  # Recompute baseline and saturated:
  
  # Form baseline model:
  if (baseline_saturated){
    x@baseline_saturated$baseline <- varcov(data,
                                            type = "chol",
                                            lowertri = "empty",
                                            vars = attr(data, "vars"),
                                            groups = attr(data, "groups"),
                                            missing = attr(data, "missing"),
                                            equal = x@equal,
                                            estimator = x@estimator,
                                            baseline_saturated = FALSE)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    x@baseline_saturated$saturated <- varcov(data,
                                             type = "chol", 
                                             lowertri = "full", 
                                             vars = attr(data, "vars"),
                                             groups = attr(data, "groups"),
                                             missing = attr(data, "missing"),
                                             equal = x@equal,
                                             estimator = x@estimator,
                                             baseline_saturated = FALSE)
    
    # if not FIML, Treat as computed:
    if (x@estimator != "FIML"){
      x@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      x@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$saturated),x@baseline_saturated$saturated)      
    }
  } else {
    x@baseline_saturated <- list()
  }
  
  
  # Return:
  return(x)
}