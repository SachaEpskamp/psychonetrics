# Inner function to make sample stats object:
samplestats <- function(
  data, # Dataset
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml"){
  
  # Check data:
  if (missing(data) & missing(covs)){
    stop("'data' and 'covs' may not both be missing")
  }
  if (!missing(data) & !missing(covs)){
    stop("'data' and 'covs' may not both *not* be missing")
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
    
    # Create covs and means arguments:
    if (nGroup == 1){
      cov <- sampleStats$cov
      class(cov) <- "matrix"
      covs <- list(as(cov,"dpoMatrix"))
      cors <- list(new("corMatrix", cov2cor(cov), sd = diag(cov)))
      means <- list(matrix(unclass(sampleStats$mean)))
    } else {
      cors <- lapply(sampleStats,function(x){
        cov <- x$cov
        class(cov) <- "matrix"
        mat <- new("corMatrix", cov2cor(cov), sd = diag(cov))
        mat
      })
      covs <- lapply(sampleStats,function(x){
        cov <- x$cov
        class(cov) <- "matrix"
        as(cov,"dpoMatrix")
      })
      means <- lapply(sampleStats,function(x){
        matrix(unclass(x$mean))
      })
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
    
    # Check if covs is array:
    if (!is.array(covs)){
      class(cors) <- "matrix"
      covs <- list(as(covs,"dpoMatrix"))
      cors <- list(new("corMatrix", cov2cor(covs), sd = diag(covs)))
    }
    # Check if means is missing:
    if (missing(means)){
      means <- lapply(1:nGroup,function(x)matrix(0,nVars,1))
    }
    
    # Check if means is matrix:
    if (!is.matrix(means)){
      means <-lapply(1:nGroup,function(x)means)  
    }
  }
  
  # Check if cov is dpoMatrix:
  for (i in seq_along(covs)){
    if (!is(covs[[i]],"dpoMatrix")){
      covs[[i]] <- as(covs[[i]], "dpoMatrix")
    }
  }
  
  # Set names:
  names(covs) <- groupNames
  names(means) <- groupNames
  
  # Generate samplestats object:
  object <- generate_psychonetrics_samplestats(covs = covs, cors = cors, means = means)
  
  # Fill groups:
  object@groups <- data.frame(
    label = groupNames,
    id = seq_along(groupNames),
    nobs = nobs, stringsAsFactors = FALSE
  )
  
  # Fill variables:
  object@variables <- data.frame(
    label = vars,
    id = seq_along(vars), stringsAsFactors = FALSE
  )
  
  # Return object:
  return(object)
}