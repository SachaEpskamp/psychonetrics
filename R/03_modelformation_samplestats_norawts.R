# Inner function to make sample stats object:
samplestats_norawts <- function(
  data, # Dataset
  vars, # character indicating the variables Extracted if missing from data - group variable
  ordered = character(0),
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = c("listwise","pairwise"),
  fimldata = FALSE,
  verbose = TRUE,
  storedata = FALSE,
  weightsmatrix = "none", #c("none","identity","full","diag")
  meanstructure = TRUE,
  covtype = c("choose","ML","UB"),
  corinput,
  standardize = c("none","z","quantile"),
  fullFIML = FALSE # <- if estimator is FIML, do full FIML. Usually not needed...
){
  missing <- match.arg(missing)
  covtype <- match.arg(covtype)
  standardize <- match.arg(standardize)
  
  # For corinput, set covtype to ML:
  if (!missing(corinput)){
    if (isTRUE(corinput)){
      if (covtype == "UB"){
        warning("Setting covtype = 'ML' because corinput = TRUE")
      }
      covtype <- "ML"
    }
  }
  # weightsmatrix <- match.arg(weightsmatrix)
  
  # Check data:
  if (missing(data) & length(ordered) >0){
    stop("Ordinal data only supported with raw data as input.")
  }
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
    
    # Standardize the data:
    if (standardize == "z"){
      for (v in seq_along(vars)){
        data[,vars[v]] <- (data[,vars[v]] - mean(data[,vars[v]],na.rm=TRUE)) / sd(data[,vars[v]],na.rm=TRUE)
      }
    } else if (standardize == "quantile"){
      for (v in seq_along(vars)){
        data[,vars[v]] <- quantiletransform(data[,vars[v]])
      }
    }
    
    # Estimate covariances:
    # lavOut <- lavaan::lavCor(data[,c(vars,groups)], missing = missing,output = "lavaan", group = groups, 
    #                          meanstructure = TRUE)
    # sampleStats <- lavaan::lavInspect(lavOut, what = "sample")
    
    # If missing is listwise, just cut out all NA rows:
    if (missing == "listwise" && !fimldata){
      data <- data[rowSums(is.na(data[,c(vars)])) == 0,]
    }
    
    # If ordered is TRUE, set all to ordered:
    if (isTRUE(ordered)){
      ordered <- vars
    }
    
    # Check if none or all are ordered:
    # FIXME: ADD POLYSERIALS LATER!!!
    
    if (!sum(vars %in% ordered) == 0 && !sum(vars %in% ordered) == length(vars)){
      stop("Either all variables or no variables may be ordered...")
    }
    # Do I need a WLS.W?
    if (length(ordered) > 0 && is.character(weightsmatrix)){
      needWLSV <- TRUE
      weightsmatrix <- list()
    } else {
      needWLSV <- FALSE
    }
    
    # Create covs and means arguments:
    if (nGroup == 1){
      if (length(ordered)  > 0){
        
        
        # Run the Cpp function:
        prepRes <- covPrepare_cpp(
          as.matrix(data[,vars]),
          vars %in% ordered,
          WLSweights = needWLSV,
          verbose = verbose
        )
        
        if (any(eigen(prepRes$covmat)$values < 0)){
          stop("Correlation matrix is not positive semi-definite.")
        }
        
        # Obtain results:
        covs <- list(as(prepRes$covmat, "matrix"))
        cors <- list(as(cov2cor(prepRes$covmat), "matrix"))
        thresholds <- list(prepRes$means_thresholds)
        means <- list(rep(NA,length(vars)))
        for (i in seq_along(vars)){
          if (!vars[i] %in% ordered){
            means[[1]][i] <- prepRes$means_thresholds[[i]]
            thresholds[[1]][[i]] <- numeric(0)
          }
        }
        
        if (needWLSV){
          weightsmatrix[[1]] <- prepRes$WLS_V
        }
        squares <- list(as(matrix(NA,nrow(prepRes$covmat),ncol(prepRes$covmat)), "matrix"))
        
        
      } else {
        cov <- (nrow(data[,c(vars)])-1)/(nrow(data[,c(vars)])) * cov(data[,c(vars)], use = switch(
          missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
        ))
        cov <- 0.5*(cov + t(cov))
        
        if (any(is.na(cov))){
         warning("NA sample covariances found. This will likely lead to erroneous estimates.")
        }
        
        # covs <- list(as(cov,"dsyMatrix"))
        covs <- list(as(cov,"matrix"))
        if (!any(is.na(cov))){
          cors <- list(as(new("corMatrix", cov2cor(cov), sd = diag(cov)),"matrix"))         
        } else {
          cors <- list()
        }
        
        dataMat <- as.matrix(data[,c(vars)])
        squares <- list(as(t(dataMat) %*% dataMat,"matrix"))
        rm(dataMat)
        
        # cors <- list(new("corMatrix", cov2cor(cov), sd = diag(cov)))
        means <- list(colMeans(data[,c(vars)], na.rm = TRUE))
        
        thresholds <- list(list())
      }
      
      # groupNames <- unique(data[[groups]])
      
    } else {
      covs <- list()
      cors <- list()
      means <- list()
      squares <- list()
      thresholds <- lapply(1:nGroup,function(x)list())
      # groupNames <- unique(data[[groups]])
      
      for (g in 1:nGroup){
        
        if (length(ordered)  > 0){
          
          # Run the Cpp function:
          prepRes <- covPrepare_cpp(
            as.matrix(data[,vars]),
            vars %in% ordered,
            WLSweights = needWLSV
          )
          
          if (any(eigen(prepRes$covmat)$values < 0)){
            stop("Correlation matrix is not positive semi-definite.")
          }
          
          # Obtain results:
          covs[[g]] <- as(prepRes$covmat, "matrix")
          cors[[g]] <- as(cov2cor(prepRes$covmat), "matrix")
          thresholds[[g]] <- prepRes$means_thresholds
          means[[g]] <- rep(NA,length(vars))
          for (i in seq_along(vars)){
            if (!vars[i] %in% ordered){
              means[[g]][i] <- prepRes$means_thresholds[[i]]
              thresholds[[g]][[i]] <- numeric(0)
            } 
          }
          
          if (needWLSV){
            weightsmatrix[[g]] <- prepRes$WLS_V
          }
          
        } else {
          
          
          
          
          subData <- data[data[[groups]] == g,c(vars)]
          cov <-  (nrow(subData)-1)/(nrow(subData)) * 
            cov(subData, use = switch(
              missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
            ))
          cov <- 0.5*(cov + t(cov))
          # covs[[g]] <- as(cov,"dsyMatrix")
          covs[[g]] <- as(cov,"matrix")
          
          if (!any(is.na(cov))){

            cors[[g]] <- as(new("corMatrix", cov2cor(cov), sd = diag(cov)),"matrix")         
          } else {
            cors[[g]] <- NA
          }
 
          subData <- as.matrix(subData)
          squares[[g]] <- as(t(subData) %*% subData,"matrix")

          means[[g]] <- colMeans(subData, na.rm = TRUE)
          
        }
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
    thresholds <- list()
    
    if (standardize != "none") warning("'standardize' ignored when raw data is not used.")
    
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
   
        
        
        
        # class(covs) <- "matrix"
        # covs <- list(as(covs,"dpoMatrix"))    
        covs <- list(as(covs,"matrix"))    
        # cors <- list(new("corMatrix", cov2cor(covs), sd = diag(covs)))
        
        if (!any(is.na(covs))){
          cors <- list(as(new("corMatrix", cov2cor(covs), sd = diag(covs)),"matrix"))
        } else {
          cors <- list(NA)
        }
        
      } else {
        cors <- lapply(covs,function(x){
          
          if (!any(is.na(x))){
            return(as(new("corMatrix", cov2cor(x), sd = diag(x)),"matrix"))
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
        covs[[g]] <- as(covsArray[,,g],"matrix")
        cors[[g]] <- as(new("corMatrix", cov2cor(covsArray[,,g]), sd = diag(covsArray[,,g])), "matrix")
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
    
    # Determine squares and covtype:
    # Maximum likelihood:
    if (covtype %in% c("ML","choose")){
      MLsquares <- list()
      for (i in seq_along(covs)){
        MLsquares[[i]] <- nobs[i] * (covs[[i]] + means[[i]] %*% t(means[[i]]))
      }
    }
    
    # Unbiased:
    if (covtype %in% c("UB","choose")){
      UBsquares <- list()
      for (i in seq_along(covs)){
        UBsquares[[i]] <- nobs[i] * (covUBtoML(covs[[i]],nobs[i]) + means[[i]] %*% t(means[[i]]))
      }
    }
    
    # If we need to choose, pick the one that looks the most like integers:
    if (covtype == "choose"){
      MLrest <- mean(round(unlist(lapply(MLsquares,as.matrix)),10) %% 1)
      UBrest <- mean(round(unlist(lapply(UBsquares,as.matrix)),10) %% 1)
      
      if (MLrest < UBrest){
        # if (verbose){
          message("Assuming denominator n was used in covariance computation (covtype = 'ML').")
        # }
        squares <- MLsquares
      } else {
        # if (verbose){
          message("Assuming denominator n-1 was used in covariance computation (covtype = 'UB').")
        # }
        squares <- UBsquares
      }
    } else if (covtype == "ML"){
      squares <- MLsquares
    } else {
      squares <- UBsquares
    }
    
    
    # Transform covs if needed:
    if (covtype == "UB"){
      for (i in seq_along(covs)){
        covs[[i]] <- covUBtoML(as.matrix(covs[[i]]), nobs[i])
      }
    }
  }
  
  # Check if cov is dpoMatrix:
  # for (i in seq_along(covs)){
  #   if (!is(covs[[i]],"dsyMatrix")){
  #     covs[[i]] <- as(covs[[i]], "dsyMatrix")
  #   }
  # }
  

  
 
  
  # Set names:
  names(covs) <- groupNames
  names(means) <- groupNames
  
  # Determine corinput (will also detect if standardized data was used as input):
  if (missing(corinput)){
    
    if (all(
      sapply(covs,function(x){
        all(abs(diag(x) - 1) < sqrt(.Machine$double.eps))
      })
    )){
      corinput <- TRUE
    } else {
      corinput <- FALSE
    }
  }
  
  # Generate samplestats object:
  object <- generate_psychonetrics_samplestats(covs = covs, cors = cors, means = means, corinput = corinput, thresholds = thresholds, squares = squares, fullFIML=fullFIML, groupvar=groups)
  
  # Fill groups:
  object@groups <- data.frame(
    label = groupNames,
    id = seq_along(groupNames),
    nobs = nobs, stringsAsFactors = FALSE
  )
  
  # Fill variables:
  object@variables <- data.frame(
    label = vars,
    id = seq_along(vars), 
    ordered = vars %in% ordered,
    
    stringsAsFactors = FALSE
  )
  
  # add fiml data (still summary statistics...):
  if (fimldata){
    if (!missing(data)){
      object@fimldata <- lapply(seq_along(groupNames),function(x){
        if (fullFIML){
          fullfimldata(data[data[[groups]] == x,vars],verbose=verbose)
        } else {
          missingpatterns(data[data[[groups]] == x,vars],verbose=verbose)  
        }
        
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
  
  # add WLS.W:
  
  if (is.list(weightsmatrix) || is.matrix(weightsmatrix)){
    if (is.list(weightsmatrix)){
      object@WLS.W <- lapply(weightsmatrix, function(x)x) # as(weightsmatrix,"Matrix")
    } else {
      object@WLS.W <- lapply(1:nGroup,function(x)weightsmatrix)
    }
    
    # Check if mean structure is added, otherwise add identitiy matrix:
    
    # FIXME: DISABLED FOR NOW
    #   if (ncol(object@WLS.W[[1]]) != nVars + nVars*(nVars+1)/2){
    #     if (ncol(object@WLS.W[[1]]) == nVars*(nVars+1)/2){
    #       if (verbose && meanstructure){
    #         warning("WLS.W only supplied for variance/covariance structure. Adding identitiy matrix to means part.")
    #       }
    #       object@WLS.W[[1]] <- rbind(
    #         cbind(Diagonal(nVars), Matrix(0,nVars,nVars*(nVars+1)/2)),
    #         cbind(Matrix(0,nVars*(nVars+1)/2,nVars), object@WLS.W[[1]])
    #         )
    #     } else {
    #       stop("WLS.W not of appropriate dimension.")  
    #     }
    #   }
    # }
  } else {
    for (g in 1:nGroup){
      if (is.character(weightsmatrix) && weightsmatrix != "none"){
        for (g in 1:nGroup){
          if (weightsmatrix == "identity"){
            object@WLS.W[[g]] <- diag(nVars + nVars*(nVars+1)/2)
          } else if (weightsmatrix == "full"){
            subData <- data[data[[groups]] == g,c(vars)]
            object@WLS.W[[g]] <- as.matrix(LS_weightsmat(subData,meanstructure=meanstructure,corinput=corinput))
          } else if (weightsmatrix == "diag"){
            subData <- data[data[[groups]] == g,c(vars)]
            object@WLS.W[[g]] <- as.matrix(LS_weightsmat(subData, type = "diagonal",meanstructure=meanstructure,corinput=corinput))
          }
        }
      }
    }
  }
  
  
  
  
  
  
  # Return object:
  return(object)
}