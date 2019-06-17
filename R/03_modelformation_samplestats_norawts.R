# Inner function to make sample stats object:
samplestats_norawts <- function(
  data, # Dataset
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = c("listwise","pairwise"),
  fimldata = FALSE,
  verbose = TRUE,
  storedata = FALSE,
  weightsmatrix = "none" #c("none","identity","full","diag")
){
  missing <- match.arg(missing)
  # weightsmatrix <- match.arg(weightsmatrix)
  
  # Check data:
  if (missing(data) & missing(covs)){
    stop("'data' and 'covs' may not both be missing")
  }
  if (!missing(data) & !missing(covs)){
    stop("'data' and 'covs' may not both *not* be missing")
  }
  if (missing(data) & storedata){
    stop("'data' may not be missing if 'storedata = TRUE'")
  }
  if (missing(data) & is.character(weightsmatrix) && weightsmatrix %in% c("full","diag")){
    stop("'data' may not be missing if estimator = 'WLS' or esitmator = 'DWLS'")
  }
  # If data is supplied:
  if (!missing(data) && !is.null(data)){
    if (!is.data.frame(data) & !is.matrix(data)){
      stop("'data' must be a data frame or matrix")
    }
    if (is.matrix(data)){
      data <- as.data.frame(data)
    }
    
    
    # If group variable is missing, add (dummy):
    if (missing(groups)|| is.null(groups)){
      groups <- "group"
      data[[groups]] <- "fullsample"
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
    
    
    # Remove all rows with full missing:
    data <- data[rowSums(is.na(data[,vars])) < nVars,]
    
    # Estimate covariances:
    # lavOut <- lavaan::lavCor(data[,c(vars,groups)], missing = missing,output = "lavaan", group = groups, 
    #                          meanstructure = TRUE)
    # sampleStats <- lavaan::lavInspect(lavOut, what = "sample")
    
    # If missing is listwise, just cut out all NA rows:
    if (missing == "listwise" && !fimldata){
      data <- data[rowSums(is.na(data[,c(vars)])) == 0,]
    }
    
    # Create covs and means arguments:
    if (nGroup == 1){
      cov <- (nrow(data[,c(vars)])-1)/(nrow(data[,c(vars)])) * cov(data[,c(vars)], use = switch(
        missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
      ))
      cov <- 0.5*(cov + t(cov))
      covs <- list(as(cov,"dsyMatrix"))
      if (!any(is.na(cov))){
        cors <- list(new("corMatrix", cov2cor(cov), sd = diag(cov)))         
      } else {
        cors <- list()
      }
      # cors <- list(new("corMatrix", cov2cor(cov), sd = diag(cov)))
      means <- list(colMeans(data[,c(vars)], na.rm = TRUE))
      # groupNames <- unique(data[[groups]])
      
    } else {
      covs <- list()
      cors <- list()
      means <- list()
      # groupNames <- unique(data[[groups]])
      
      for (g in 1:nGroup){
        subData <- data[data[[groups]] == g,c(vars)]
        cov <-  (nrow(subData)-1)/(nrow(subData)) * 
          cov(subData, use = switch(
            missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
          ))
        cov <- 0.5*(cov + t(cov))
        covs[[g]] <- as(cov,"dsyMatrix")
        
        if (!any(is.na(cov))){
          cors[[g]] <- new("corMatrix", cov2cor(cov), sd = diag(cov))          
        } else {
          cors[[g]] <- NA
        }
        
        means[[g]] <- colMeans(subData, na.rm = TRUE)
      }
      
      
      # cors <- lapply(sampleStats,function(x){
      #   cov <- x$cov
      #   class(cov) <- "matrix"
      #   mat <- new("corMatrix", cov2cor(cov), sd = diag(cov))
      #   mat
      # })
      # covs <- lapply(sampleStats,function(x){
      #   cov <- x$cov
      #   class(cov) <- "matrix"
      #   as(cov,"dpoMatrix")
      # })
      # means <- lapply(sampleStats,function(x){
      #   matrix(unclass(x$mean))
      # })
      # groupNames <- names(sampleStats)
    }
    if (!missing(nobs)){
      warning("'nobs' argument ignored and obtained from data")
    }
    
    nobs <- as.vector(tapply(data[[groups]],data[[groups]],length))
  } else {
    ### Input via matrices ###
    # Check groups:
    if (missing(groups) || is.null(groups)){
      if (is.array(covs) && length(dim(covs)) > 2){
        if (!is.null(dimnames(covs)[[3]])){
          groups <- groupNames <- paste0("group_",seq_len(dim(covs)[[3]]))  
        } else {
          groups <- groupNames <- dimnames(covs)[[3]]
        }
      } else if (is.list(covs)){
        if (!is.null(names(covs))){
          groups <- groupNames <- names(covs)
        } else {
          groups <- groupNames <- paste0("group_",seq_len(length(covs)))  
        }
      } else {
        groups <- groupNames <- "fullsample" 
      }
    } else {
      groupNames <- groups
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
      if (!is.list(covs)){
        class(covs) <- "matrix"
        covs <- list(as(covs,"dpoMatrix"))        
        # cors <- list(new("corMatrix", cov2cor(covs), sd = diag(covs)))
        
        if (!any(is.na(covs))){
          cors <- list(new("corMatrix", cov2cor(covs), sd = diag(covs)))
        } else {
          cors <- list(NA)
        }
        
      } else {
        cors <- lapply(covs,function(x){
          
          if (!any(is.na(x))){
            return(new("corMatrix", cov2cor(x), sd = diag(x)))
          } else {
            return(NA)
          }
          
          # new("corMatrix", cov2cor(x), sd = diag(x))
        })
      }
      
      nVars <- ncol(covs[[1]])
    } else {
      
      # Make array
      if (length(dim(covs)) == 2){
        if (!is.null(colnames(covs))){
          vars <- colnames(covs)
        } else {
          vars <- paste0("V",seq_len(ncol(covs)))
        }
        covs <- array(covs,c(dim(covs),nGroup))
        dimnames(covs) <- list(vars,vars,NULL)
      }
      
      # Now create list:
      covsArray <- covs
      covs <- list()
      cors <- list()
      for (g in 1:nGroup){
        covs[[g]] <- as(covsArray[,,g],"dpoMatrix")
        cors[[g]] <- new("corMatrix", cov2cor(covsArray[,,g]), sd = diag(covsArray[,,g])) 
      }
    }
    
    # Number of vars:
    nVars <- nrow(covs[[1]])
    
    # Varnames:
    if (missing(vars)){
      if (!is.null(colnames(covs[[1]]))){
        vars <- colnames(covs[[1]])
      } else {
        vars <- paste0("V",seq_len(nVars))
      }
    }
    
    # Check if means is missing:
    if (missing(means)){
      means <- lapply(1:nGroup,function(x)matrix(0,nVars,1))
    }
    
    # Check if means is matrix:
    if (!is.list(means)){
      means <-lapply(1:nGroup,function(x)means)  
    }
  }
  
  # Check if cov is dpoMatrix:
  for (i in seq_along(covs)){
    if (!is(covs[[i]],"dsyMatrix")){
      covs[[i]] <- as(covs[[i]], "dsyMatrix")
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
  
  # add fiml data (still summary statistics...):
  if (fimldata){
    if (!missing(data)){
      object@fimldata <- lapply(seq_along(groupNames),function(x){
        missingpatterns(data[data[[groups]] == x,vars],verbose=verbose)
      })
    } else {
      object@fimldata <- lapply(seq_along(groupNames),function(x){
        missingpatterns_covs(means[[x]],covs[[x]],nobs[x],verbose=verbose)
      })
    }
 
  }
  
  # add full data:
  if (storedata){
    # Overwrite group with name:
    data[[groups]] <- groupNames[ data[[groups]]]
    
    object@rawdata <- data[,c(vars, groups)]
    attr(object@rawdata, "vars") <- vars
    attr(object@rawdata, "groups") <- groups
    attr(object@rawdata, "missing") <- missing
    attr(object@rawdata, "fimldata") <- fimldata
  }
  
  # add WLS.V:
  if (is.list(weightsmatrix) || is.matrix(weightsmatrix)){
    if (is.list(weightsmatrix)){
      object@WLS.V <- as(weightsmatrix,"Matrix")
    } else {
      object@WLS.V <- lapply(1:nGroup,function(x)as(weightsmatrix,"Matrix"))
    }
    
    # Check if mean structure is added, otherwise add identitiy matrix:
    for (g in 1:nGroup){
      if (ncol(object@WLS.V[[1]]) != nVars + nVars*(nVars+1)/2){
        if (ncol(object@WLS.V[[1]]) == nVars*(nVars+1)/2){
          if (verbose){
            warning("WLS.V only supplied for variance/covariance structure. Adding identitiy matrix to means part.")
          }
          object@WLS.V[[1]] <- rbind(
            cbind(Diagonal(nVars), Matrix(0,nVars,nVars*(nVars+1)/2)),
            cbind(Matrix(0,nVars*(nVars+1)/2,nVars), object@WLS.V[[1]])
            )
        } else {
          stop("WLS.V not of appropriate dimension.")  
        }
      }
    }
  } else {
    if (weightsmatrix != "none"){
      for (g in 1:nGroup){
        if (weightsmatrix == "identity"){
          object@WLS.V[[g]] <- Diagonal(nVars + nVars*(nVars+1)/2)  
        } else if (weightsmatrix == "full"){
          subData <- data[data[[groups]] == g,c(vars)]
          object@WLS.V[[g]] <- LS_weightsmat(subData)
        } else if (weightsmatrix == "diag"){
          subData <- data[data[[groups]] == g,c(vars)]
          object@WLS.V[[g]] <- LS_weightsmat(subData, type = "diagonal")
        }
      }
    }
  }
  
  
  
  
  # Return object:
  return(object)
}