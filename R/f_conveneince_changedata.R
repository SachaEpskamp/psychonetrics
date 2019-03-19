changedata <- function(x, data, covs, nobs, means, groups, missing = "listwise"){
  if (x@distribution != "Gaussian"){
    stop("Only Gaussian distirbution models supported.")
  }
  if (x@model == "var1"){
    stop("VAR models not supported yet.")
  }
  
  
 # Set computed to false:
  x@computed <- FALSE
  
  # Check if labels are the same:
  if (!all(x@sample@variables$label %in% colnames(data))){
    stop("New dataset does not contain the same variables as used in model.")
  }
  
  # Variables:
  vars <-  x@sample@variables$label
  
  # Estimator:
  estimator <- x@estimator
  
  # Data set:
  x@sample <- samplestats(data = data, 
              vars = vars, 
              groups = groups,
              covs = covs, 
              means = means, 
              nobs = nobs, 
              missing  = ifelse(estimator == "FIML","pairwise",missing),
              fimldata = estimator == "FIML")
  
  # Number of observations:
  nVar <- nrow(x@sample@variables)
  nGroup <- nrow(x@sample@groups)
  x@sample@nobs <-  
    nVar * (nVar+1) / 2 * nGroup + # Covariances per group
    nVar * nGroup
  
  # Add baseline and saturated:
  
  # Form baseline model:
  x@baseline_saturated$baseline <- varcov(data,
                                              type = "chol",
                                              lowertri = "empty",
                                              vars = vars,
                                              groups = groups,
                                              covs = covs,
                                              means = means,
                                              nobs = nobs,
                                              missing = missing,
                                              equal = x@equal,
                                              estimator = estimator,
                                              baseline_saturated = FALSE)
  
  # Add model:
  # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
  
  
  ### Saturated model ###
  x@baseline_saturated$saturated <- varcov(data,
                                               type = "chol", 
                                               lowertri = "full", 
                                               vars = vars,
                                               groups = groups,
                                               covs = covs,
                                               means = means,
                                               nobs = nobs,
                                               missing = missing,
                                               equal = x@equal,
                                               estimator = estimator,
                                               baseline_saturated = FALSE)
  
  # if not FIML, Treat as computed:
  if (estimator != "FIML"){
    x@baseline_saturated$saturated@computed <- TRUE
    
    # FIXME: TODO
    x@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$saturated),x@baseline_saturated$saturated)      
  }
  
  return(x)
}
