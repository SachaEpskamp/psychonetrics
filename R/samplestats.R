# Inner function to make sample stats object:
samplestats <- function(
  data, # Dataset
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covmat, # alternative covmat (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml"){
  
  # Check data:
  if (missing(data) & missing(covmat)){
    stop("'data' and 'covmat' may not both be missing")
  }
  if (!missing(data) & !missing(covmat)){
    stop("'data' and 'covmat' may not both *not* be missing")
  }
  
  # If data is supplied:
  if (!missing(data)){
    if (!is.data.frame(data) & !is.matrix(data)){
      stop("'data' must be a data frame or matrix")
    }
    if (is.matrix(data)){
      data <- as.data.frame(data)
    }
    # If group variable is missing, add (dummy):
    if (missing(groups)){
      groups <- "singlegroup"
      data[[groups]] <- "singlegroup"
    }
    
    # Extract group names:
    groupNames <- unique(data[[groups]])
    
    # number of groups:
    nGroup <- length(groupNames)
    
    # Overwrite group with integer:
    data[[groups]] <- match(data[[groups]], groupNames)
    
    # If vars is missing, obtain from data:
    if (missing(vars)){
      vars <- names(data[,names(data)!=groups])
    }
    
    # Number of variables:
    nVars <- length(vars)
    
    # Estimate covariances:
    lavOut <- lavaan::lavCor(data[,c(vars,groups)], missing = missing,output = "lavaan", group = groups)
    sampleStats <- lavaan::lavInspect(lavOut, what = "sample")
    
    # Create covmat and means arguments:
    if (nGroup == 1){
      covmat <- array(sampleStats$cov, c(nVars,nVars,1))
      means <- matrix(sampleStats$mean, nVars, 1)
    } else {
      covmat <- do.call(abind::abind, c(lapply(sampleStats,"[[","cov"), list(along = 3)))
      means <- do.call(abind::abind, c(lapply(sampleStats,"[[","mean"), list(along = 2)))
      groupNames <- names(sampleStats)
    }
    if (!missing(nobs)){
      warning("'nobs' argument ignored and obtained from data")
    }
    nobs <- lavaan::lavInspect(lavOut, "nobs")
  } else {
    ### Input via matrices ###
    # Check groups:
    if (missing(groups)){
      groups <- groupNames <- "singlegroup"
    }
    nGroup <- length(groups)
    
    # if nobs missing, stop:
    if (missing(nobs)){
      stop("'nobs' may not be missing")
    }
    if (length(nobs) != nGroup){
      stop("'nobs' must be a vector with sample size per group")
    }
    
    # Check if covmat is array:
    if (!is.array(covmat)){
      covmat <- array(covmat, c(nVars, nVars, nGroup))
    }
    # Check if means is missing:
    if (missing(means)){
      means <- matrix(0, nVars, nGroup)
    }
    
    # Check if means is matrix:
    if (!is.matrix(means)){
      means <- matrix(means, nVars, nGroup)
    }
  }
  
  # Generate samplestats object:
  object <- generate_psychonetrics_samplestats(covmat = covmat,means = means)
  
  # Fill groups:
  object@groups <- data.frame(
    label = groupNames,
    id = seq_along(groupNames),
    nobs = nobs
  )
  
  # Fill variables:
  object@variables <- data.frame(
    label = vars,
    id = seq_along(vars)
  )
  
  # Return object:
  return(object)
}