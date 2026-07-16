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
  fullFIML = FALSE, # <- if estimator is FIML, do full FIML. Usually not needed...
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  likelihood = c("normal","wishart"), # Gaussian likelihood scaling. "wishart" stores the UNBIASED (n-1 denominator) sample covariance per group instead of the ML (n denominator) one (matches lavaan likelihood = "wishart").
  compute_ml_gamma = FALSE, # Compute + store the asymptotic covariance Gamma of the sample statistics for a complete-data ML fit. Only needed by the robust.sem estimators (MLM/MLMV/MLMVS); off by default so plain ML fits do not pay the (potentially large) fourth-order-moment cost. Set TRUE by the constructors when a robust.sem estimator is requested.
  weights = NULL # Optional numeric vector of SAMPLING weights (pseudo-ML), aligned with the rows of 'data'. When supplied, each group's moments are computed as WEIGHTED means/covariances (weights normalized within group to sum to n_g) and the normalized per-group weights are stored for the robust (MLR) sandwich standard errors. Only for complete-data raw input.
){
  likelihood <- match.arg(likelihood)
  # Per-group normalized sampling weights (filled below when 'weights' is given):
  sampling_weights_list <- list()
  # bootstrap defaults:
  if (bootstrap == "nonparametric"){
    bootstrap <- TRUE
    
    if (!missing(boot_sub)){
       warning("bootstrap = 'nonparametric' overwrites boot_sub = 1")
    }
    boot_sub <- 1
    
    if (!missing(boot_resample)){
      warning("bootstrap = 'nonparametric' overwrites boot_resample = TRUE")
    }
    boot_resample <- TRUE
  } else if (bootstrap == "case"){
    bootstrap <- TRUE
    
    if (!missing(boot_sub)){
      warning("bootstrap = 'case' overwrites boot_sub = 0.75")
    }
    boot_sub <- 0.75
    
    if (!missing(boot_resample)){
      warning("bootstrap = 'case' overwrites boot_resample = FALSE")
    }
    boot_resample <- FALSE
  } else if (!is.logical(bootstrap)){
    stop("'bootstrap' must be TRUE, FALSE, 'nonparametric' or 'case'")
  } else {
    if (missing(boot_sub)){
      boot_sub <- 1
    }
    if (missing(boot_resample)){
      boot_resample <- TRUE
    }
    
    if (boot_sub < 0 || boot_sub > 1){
      stop("'boot_sub' must be between 0 and 1.")
    }
    
    if (!is.logical(boot_resample)){
      stop("'boot_resample' must be logical.")
    }
  }
  
  
  # Error if missing data for bootstrap:
  if (isTRUE(bootstrap) && missing(data)){
    stop("Bootstrap requires data to be supplied")
  }
    
    
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

  # Initialize Gamma matrix storage for WLSMV (must be before data/covs branching):
  gammamatrix <- list()

  # If data is supplied:
  if (!missing(data) && !is.null(data)){
    if (!is.data.frame(data) & !is.matrix(data)){
      stop("'data' must be a data frame or matrix")
    }
    if (is.matrix(data)){
      data <- as.data.frame(data)
    }

    # Ensure plain data.frame (tibbles have different subsetting behavior):
    if (inherits(data, "tbl_df")){
      data <- as.data.frame(data)
    }


    # If group variable is missing, add (dummy):
    if (missing(groups)|| is.null(groups)){
      groups <- "group"
      data[[groups]] <- "fullsample"
    } else {
      if (!groups %in% names(data)){
        stop("'groups' object does not refer to a column in the data")
      }
      
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
    
    
    # Remove all rows with full missing (drop = FALSE keeps a single-variable
    # selection a data frame so rowSums() works; e.g. a bivariate meta_varcov
    # has only one meta-level variable):
    data <- data[rowSums(is.na(data[,vars,drop=FALSE])) < nVars,]
    
    ### Bootstrap the data ###
    
    if (isTRUE(bootstrap)){
      
      data <- data[sample(seq_len(nrow(data)), round(boot_sub*nrow(data)), replace = boot_resample),]
      
    }
    
    ###
    
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
      data <- data[rowSums(is.na(data[,c(vars),drop=FALSE])) == 0,]
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
      ordinal_wtype <- weightsmatrix  # Preserve: "full", "diag", or "identity"
      weightsmatrix <- list()
    } else {
      needWLSV <- FALSE
      ordinal_wtype <- "none"
    }

    # Early, informative checks for the continuous (cov()) path. Variables that
    # are not declared 'ordered' must be numeric, otherwise cov() fails with a
    # cryptic "x must be numeric" error (e.g. when factor/character Ising data is
    # accidentally passed to a continuous model):
    contVars <- vars[!(vars %in% ordered)]
    if (length(contVars) > 0){
      nonNumeric <- contVars[!vapply(contVars, function(v) is.numeric(data[[v]]), logical(1))]
      if (length(nonNumeric) > 0){
        stop(paste0("The following variable(s) are not numeric: ", paste0(nonNumeric, collapse = ", "),
                    ". Non-numeric (e.g. factor or character) variables must be declared in 'ordered', or converted to numeric (for example for Ising models, recode responses to integers)."))
      }
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
          # WLS_V from covPrepare_cpp is full Gamma inverse
          # For full WLS: store as-is (full weight matrix)
          # For DWLS: use 1/diag(Gamma), NOT diag(inverse(Gamma)) -- matches lavaan convention
          Gamma_full <- solve(as.matrix(prepRes$WLS_V))
          gammamatrix[[1]] <- Gamma_full

          if (ordinal_wtype == "diag"){
            # DWLS: diagonal weight = 1/diag(Gamma)
            weightsmatrix[[1]] <- diag(1 / diag(Gamma_full))
          } else {
            # WLS (full): use the full inverse of Gamma
            weightsmatrix[[1]] <- prepRes$WLS_V
          }
        }
        squares <- list(as(matrix(NA,nrow(prepRes$covmat),ncol(prepRes$covmat)), "matrix"))
        
        
      } else {
        # drop = FALSE keeps a single-variable selection a matrix so nrow()/cov()
        # behave (otherwise a 1-node model fails with a cryptic NULL nrow error):
        datmat <- data[,c(vars), drop = FALSE]
        # ML (n denominator) covariance by default; the UNBIASED (n-1) covariance
        # under likelihood = "wishart" (cov() is already n-1, so the rescaling is
        # skipped):
        cov_rescale <- if (likelihood == "wishart") 1 else (nrow(datmat)-1)/(nrow(datmat))
        cov <- cov_rescale * cov(datmat, use = switch(
          missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
        ))

        # For n=1, make the covariances 0:
        if (nrow(data)==1){
          cov[is.na(cov)] <- 0
        }
        
        cov <- 0.5*(cov + t(cov))
        
        if (any(is.na(cov))){
         warning("NA sample covariances found. This will likely lead to erroneous estimates.")
        }
        
        # covs <- list(as(cov,"dsyMatrix"))
        covs <- list(as(cov,"matrix"))
        if (!any(is.na(cov))){
          cors <- list(as(as(cov2cor(cov), "symmetricMatrix"),"matrix"))         
        } else {
          cors <- list()
        }
        
        dataMat <- as.matrix(data[,c(vars),drop=FALSE])
        squares <- list(as(t(dataMat) %*% dataMat,"matrix"))
        rm(dataMat)

        # cors <- list(new("corMatrix", cov2cor(cov), sd = diag(cov)))
        means <- list(colMeans(data[,c(vars),drop=FALSE], na.rm = TRUE))
        
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

          # Subset the rows belonging to this group (data[[groups]] has been
          # recoded to integer group ids 1:nGroup above), so that each group
          # gets its own polychorics/thresholds/weights instead of the pooled
          # sample:
          subData <- data[data[[groups]] == g, vars]

          # Run the Cpp function:
          prepRes <- covPrepare_cpp(
            as.matrix(subData),
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

          # Squares are not used for ordinal data; assign an NA matrix for
          # symmetry with the single-group branch:
          squares[[g]] <- as(matrix(NA,nrow(prepRes$covmat),ncol(prepRes$covmat)), "matrix")

          if (needWLSV){
            Gamma_full_g <- solve(as.matrix(prepRes$WLS_V))
            gammamatrix[[g]] <- Gamma_full_g

            if (ordinal_wtype == "diag"){
              weightsmatrix[[g]] <- diag(1 / diag(Gamma_full_g))
            } else {
              weightsmatrix[[g]] <- prepRes$WLS_V
            }
          }

        } else {
          
          
          
          
          # drop = FALSE keeps a single-variable selection a matrix (see above):
          subData <- data[data[[groups]] == g,c(vars), drop = FALSE]
          # ML (n) covariance by default; unbiased (n-1) under wishart:
          cov_rescale <- if (likelihood == "wishart") 1 else (nrow(subData)-1)/(nrow(subData))
          cov <-  cov_rescale *
            cov(subData, use = switch(
              missing, "listwise" = "complete.obs", "pairwise" = "pairwise.complete.obs"
            ))
          cov <- 0.5*(cov + t(cov))
          # covs[[g]] <- as(cov,"dsyMatrix")
          covs[[g]] <- as(cov,"matrix")
          
          if (!any(is.na(cov))){

            cors[[g]] <- as(as(cov2cor(cov), "symmetricMatrix"),"matrix")         
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

    # Sampling weights (pseudo-ML): recompute each group's moments as WEIGHTED
    # means / covariances and record the normalized per-group weights (each
    # summing to n_g, lavaan's sampling.weights.normalization = "total"). Only
    # the continuous complete-data path is supported (the constructors guard
    # against ordinal / missing / summary-statistic input).
    if (!is.null(weights)){
      weights <- as.numeric(weights)
      if (length(weights) != nrow(data)){
        stop("Internal error: sampling weights length does not match the number of data rows.")
      }
      if (any(!is.finite(weights)) || any(weights < 0)){
        stop("Sampling weights must be finite and non-negative.")
      }
      if (length(ordered) > 0){
        stop("Sampling weights are not (yet) supported for ordinal data.")
      }
      # Normalize GLOBALLY so the weights over ALL groups sum to the total sample
      # size (lavaan's default sampling.weights.normalization = "total"). The
      # moments are ratios and so are invariant to the normalization, but the
      # sandwich meat depends on the absolute weight values sum_i w_i^2 s_i s_i',
      # so the absolute scale must match lavaan for multi-group models.
      Ntot <- nrow(data)
      wtot <- sum(weights)
      if (wtot <= 0) stop("Sampling weights sum to zero.")
      weights <- weights * (Ntot / wtot)
      grpvec <- data[[groups]]
      for (g in 1:nGroup){
        idx <- which(grpvec == g)
        w <- weights[idx]
        swg <- sum(w)
        if (swg <= 0) stop("Sampling weights within a group sum to zero.")
        Y <- as.matrix(data[idx, vars, drop = FALSE])
        if (anyNA(Y)){
          stop("Sampling weights are not (yet) supported with missing data.")
        }
        mw <- as.numeric(crossprod(w, Y) / swg)            # weighted mean (divisor = sum of group weights)
        Yc <- sweep(Y, 2, mw)
        Sw <- crossprod(Yc * sqrt(w)) / swg                # weighted ML covariance
        dn <- list(vars, vars)
        dimnames(Sw) <- dn
        covs[[g]] <- as(Sw, "matrix")
        means[[g]] <- stats::setNames(mw, vars)
        cors[[g]] <- if (!any(is.na(Sw))) as(cov2cor(Sw), "matrix") else matrix(NA, length(vars), length(vars))
        squares[[g]] <- as(crossprod(Y * sqrt(w)) / swg * length(idx), "matrix")   # weighted raw sums of squares
        sampling_weights_list[[g]] <- w
        # The pseudo-log-likelihood weights each group by its total weight sum
        # (sum_i w_i), not its raw count, so set the group's effective sample
        # size to swg. This makes the fit function, gradient, expected/observed
        # information and chi-square (= sum_g swg * F_g) match lavaan; the total
        # sum_g swg = N (global normalization) keeps n = N for BIC etc.
        nobs[g] <- swg
      }
    }
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
          cors <- list(as(as(cov2cor(covs), "symmetricMatrix"),"matrix"))
        } else {
          cors <- list(NA)
        }
        
      } else {
        cors <- lapply(covs,function(x){
          
          if (!any(is.na(x))){
            return(as(as(cov2cor(x), "symmetricMatrix"),"matrix"))
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
        cors[[g]] <- as(as(cov2cor(covsArray[,,g]), "symmetricMatrix"), "matrix")
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

      # The "choose" heuristic relies on the sufficient statistics (squares)
      # being integer-valued (as they are for e.g. Ising count data). If neither
      # the ML nor the UB scaling yields near-integer squares, the input
      # responses are likely non-integer and the guess is unreliable:
      if (min(MLrest, UBrest) > sqrt(.Machine$double.eps)){
        warning("Could not reliably determine the covariance denominator (covtype) from the input: the implied sufficient statistics are not integer-valued. Please set 'covtype' explicitly ('ML' or 'UB').")
      }

      if (MLrest < UBrest){
        # if (verbose){
          message("Assuming denominator n was used in covariance computation (covtype = 'ML').")
        # }
        squares <- MLsquares
        covtype <- "ML"
      } else {
        # if (verbose){
          message("Assuming denominator n-1 was used in covariance computation (covtype = 'UB').")
        # }
        squares <- UBsquares
        covtype <- "UB"
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

    # Under likelihood = "wishart" the stored covariance must be the UNBIASED
    # (n-1 denominator) one. The covs are in ML form here, so rescale by
    # n/(n-1) per group:
    if (likelihood == "wishart"){
      for (i in seq_along(covs)){
        covs[[i]] <- as.matrix(covs[[i]]) * nobs[i] / (nobs[i] - 1)
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
  
  # Determine corinput (will also detect if standardized data was used as
  # input). isTRUE(): with heavy missingness a (pairwise) covariance can carry
  # NA diagonal elements, which must read as "not correlation input" rather
  # than crash the check (the NA-covariance warning above already flags them):
  if (missing(corinput)){
    if (isTRUE(all(
      sapply(covs,function(x){
        all(abs(diag(x) - 1) < sqrt(.Machine$double.eps))
      })
    ))){
      corinput <- TRUE
    } else {
      corinput <- FALSE
    }
  }
  
  # Generate samplestats object:
  object <- generate_psychonetrics_samplestats(covs = covs, cors = cors, means = means, corinput = corinput, 
                                               thresholds = thresholds, squares = squares, fullFIML=fullFIML, 
                                               groupvar=groups, bootstrap = bootstrap, boot_sub = boot_sub,
                                               boot_resample = boot_resample)
  
  # Store the per-group normalized sampling weights (empty list if none used):
  object@weights <- sampling_weights_list

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
          fullfimldata(data[data[[groups]] == x,vars,drop=FALSE],verbose=verbose)
        } else {
          missingpatterns(data[data[[groups]] == x,vars,drop=FALSE],verbose=verbose)
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
            nW <- meanstructure * nVars + nVars*(nVars+1)/2 - corinput * nVars
            object@WLS.W[[g]] <- diag(nW)
            # Store full Gamma so robust (sandwich) SEs are available for the
            # ULS estimator. Only possible when raw data is supplied; with a
            # covariance matrix as input Gamma cannot be computed and the SEs
            # fall back to the non-robust form (with a warning) in getVCOV().
            if (!missing(data) && !is.null(data)){
              subData <- data[data[[groups]] == g,c(vars),drop=FALSE]
              gammamatrix[[g]] <- LS_Gamma(subData, meanstructure=meanstructure, corinput=corinput)
            }
          } else if (weightsmatrix == "full"){
            subData <- data[data[[groups]] == g,c(vars),drop=FALSE]
            object@WLS.W[[g]] <- as.matrix(LS_weightsmat(subData,meanstructure=meanstructure,corinput=corinput))
            # Store full Gamma for WLSMV correction (can be obtained by inverting full W):
            gammamatrix[[g]] <- solve(object@WLS.W[[g]])
          } else if (weightsmatrix == "diag"){
            subData <- data[data[[groups]] == g,c(vars),drop=FALSE]
            object@WLS.W[[g]] <- as.matrix(LS_weightsmat(subData, type = "diagonal",meanstructure=meanstructure,corinput=corinput))
            # Store full Gamma for WLSMV correction (cannot recover from diagonal W, compute from scratch):
            gammamatrix[[g]] <- LS_Gamma(subData, meanstructure=meanstructure, corinput=corinput)
          }
        }
      }
    }
  }
  
  
  
  
  
  
  # Store Gamma (asymptotic covariance) of the sample statistics for the
  # maximum-likelihood estimators (estimator = "ML") when raw, complete data
  # are available. This is needed for the robust ML estimators MLM/MLMV/MLMVS
  # (Browne 1984 sandwich SEs and the Satorra-Bentler family of scaled test
  # statistics), which map internally to estimator = "ML". The same LS_Gamma()
  # used for WLS/DWLS is applied here, honouring meanstructure / corinput row
  # dropping. Only computed when the default weightsmatrix = "none" path is
  # taken (i.e. not WLS/DWLS/ULS, which fill gammamatrix above) and the data are
  # complete (no missingness); FIML (missing data) robust ML is Phase 2.
  # Gated on compute_ml_gamma so plain (non-robust) ML fits do NOT pay the
  # fourth-order-moment cost (Gamma is O(p^4) in memory/time -- e.g. ~84 MB and
  # ~12 s at p = 80); the constructors request it only for MLM/MLMV/MLMVS, and
  # setestimator() computes it lazily when switching a stored-data model to a
  # robust.sem estimator.
  if (isTRUE(compute_ml_gamma) &&
      length(gammamatrix) == 0 &&
      is.character(weightsmatrix) && length(weightsmatrix) == 1 && weightsmatrix == "none" &&
      !missing(data) && !is.null(data)){
    ml_gamma <- tryCatch({
      gm <- vector("list", nGroup)
      ok <- TRUE
      grp <- data[[groups]]
      for (g in seq_len(nGroup)){
        # The group column is integer-coded here, EXCEPT when storedata = TRUE
        # has already overwritten it with the group names (see the storedata
        # block above). Match either representation:
        mask <- grp == g | grp == groupNames[g]
        subData <- data[mask, c(vars), drop = FALSE]
        # Only well-defined for complete data (listwise): skip if any NA remain.
        if (anyNA(subData) || nrow(subData) < 2){ ok <- FALSE; break }
        gm[[g]] <- LS_Gamma(subData, meanstructure = meanstructure, corinput = corinput)
      }
      if (ok) gm else list()
    }, error = function(e) list())
    if (length(ml_gamma) == nGroup){
      gammamatrix <- ml_gamma
    }
  }

  # Store Gamma (asymptotic covariance) for WLSMV correction:
  if (length(gammamatrix) > 0){
    object@WLS.Gamma <- gammamatrix
  }

  # Return object:
  return(object)
}