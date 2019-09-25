# Latent network model creator
meta_varcov <- function(
  cors, # List of correlation matrices as input. Must contain NAs
  nobs, # vector of sample sizes as input
  
  Vmats, # Optional list of V matrices for each group. Will be averaged. 
  Vmethod = c("default","psychonetrics_individual", "psychonetrics_pooled", "metaSEM_individual","metaSEM_weighted"), # How to obtain V matrices if Vmats is not supplied?
  
  # Model setup:
  type = c("cor", "ggm"), # Same as in varcov. Currently only cor and ggm are supported.
  sigma = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa = "full", # Precision
  # rho = "full", # Correlations
  omega = "full", # Partial correlations
  lowertri = "full", # Cholesky
  delta = "full", # Used for both ggm and pcor
  rho = "full", # Used for cor
  SD = "full", # Used for cor
  
  # Random effects setup:
  randomEffects = c("cov","chol","prec","ggm","cor"),
  sigma_randomEffects = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa_randomEffects = "full", # Precision
  # rho = "full", # Correlations
  omega_randomEffects = "full", # Partial correlations
  lowertri_randomEffects = "full", # Cholesky
  delta_randomEffects = "full", # Used for both ggm and pcor
  rho_randomEffects = "full", # Used for cor
  SD_randomEffects = "full", # Used for cor
  
  vars, # Only used for labeling the variables
  
  # Some extra stuff:
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  optimizer = "default",
  
  sampleStats, # Leave to missing
  verbose = TRUE
){
  # For now, I will always assume correlations were used in the input:
  corinput <- TRUE
  
  randomEffects <- match.arg(randomEffects)
  type <- match.arg(type)
  Vmethod <- match.arg(Vmethod)
  if (Vmethod == "default"){
    Vmethod <- "psychonetrics_individual"
    if (!missing(Vmats)){
      stop("'Vmats' must be missing if 'Vmethod' is not 'default'")
    }
  }
  
  # set the labels:
  if (missing(vars)){
    vars <- unique(unlist(sapply(cors,colnames)))
    if (is.null(vars)){
      vars <- paste0("V",seq_len(max(sapply(cors,ncol))))  
      if (length(unique(sapply(cors,colnames))) > 1){
        stop("Correlation matrices of different dimensions not supported without column labels.")
      }
      for (i in seq_along(cors)){
        rownames(cors[[i]]) <- colnames(cors[[i]]) <- vars
      }
    }
  } else {
    if (is.null(colnames(cors[[1]]))){
      for (i in seq_along(cors)){
        rownames(cors[[i]]) <- colnames(cors[[i]]) <- vars
      }
    }
  }
  
  # Number of nodes:
  nNode <- length(vars)
  
  # Form the dataset from the lower triangles of cor matrices:
  data <- dplyr::bind_rows(lapply(cors,function(x){
    df <- as.data.frame(t(x[lower.tri(x)]))
    names(df) <- paste0(rownames(x)[row(x)[lower.tri(x)]], " -- ",colnames(x)[col(x)[lower.tri(x)]])
    df
  }))
  
  # Labels for the correlations:
  corvars <- colnames(data)
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    
    sampleStats <- samplestats(data = data, 
                               vars = corvars, 
                               missing = "pairwise",
                               fimldata = TRUE,
                               storedata = FALSE,
                               meanstructure = TRUE,
                               verbose=verbose)
  }
  
  # Overwrite corInput:
  sampleStats@corinput <- corinput
  
  
  # Check some things:
  nCor <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "meta_varcov", sample = sampleStats, computed = FALSE, 
                                  optimizer = optimizer, estimator = "FIML", distribution = "Gaussian",
                                  types = list(y = type, randomEffects = randomEffects),
                                  submodel = type, meanstructure = TRUE)
  
  # Number of groups:
  nGroup <- 1
  
  
  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))
  
  # Add number of observations:
  model@sample@nobs <-  
    nCor * (nCor-1) / 2 * nGroup + # Covariances per group
    nCor * nGroup + # Variances (ignored if correlation matrix is input)
    nMeans 
  
  
  ### Estimate V matrices
  if (!missing(Vmats)){
    avgVmat <- Reduce("+", Vmats) / length(Vmats)
    
  } else {
    
    if (verbose){
      message("Computing sampling error approximation...")
    }
    
    if (Vmethod == "psychonetrics_individual"){
      # For each group, make a model and obtain VCOV:
      Vmats <- lapply(seq_along(cors),function(i){
        mod <- varcov(covs=cors[[i]],nobs=sampleSizes[i], corinput = TRUE, type = "cor", baseline_saturated = FALSE, verbose = FALSE)
        mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
        getVCOV(mod)
      })
      
      avgVmat <- Reduce("+", Vmats) / length(Vmats)
      
    }
    
    
    if (Vmethod == "psychonetrics_pooled"){
      # Single multi-group model:
      mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = "rho", baseline_saturated = FALSE, verbose = FALSE)
      mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
      acov <- getVCOV(mod)
      avgVmat <- acov * length(cors)
    }
    
    if (Vmethod == "metaSEM_individual"){
      acovs <- metaSEM::asyCov(cors, sampleSizes, acov = "individual")
      vMats <- list()
      for (i in seq_len(nrow(acovs))){
        vMats[[i]] <- matrix(0, nCor, nCor)
        vMats[[i]][lower.tri(vMats[[i]],diag=TRUE)] <- acovs[i,]
        vMats[[i]][upper.tri(vMats[[i]],diag=TRUE)] <- t(vMats)[[i]][upper.tri(vMats[[i]],diag=TRUE)]
      }
      avgVmat <- Reduce("+", Vmats) / length(Vmats)
    }
    
    
    
    if (Vmethod == "metaSEM_weighted"){
      acovs <- metaSEM::asyCov(cors, sampleSizes, acov = "weighted")
      vMats <- list()
      for (i in seq_len(nrow(acovs))){
        vMats[[i]] <- matrix(0, nCor, nCor)
        vMats[[i]][lower.tri(vMats[[i]],diag=TRUE)] <- acovs[i,]
        vMats[[i]][upper.tri(vMats[[i]],diag=TRUE)] <- t(vMats)[[i]][upper.tri(vMats[[i]],diag=TRUE)]
      }
      avgVmat <- Reduce("+", Vmats) / length(Vmats)
    }
  }
  
  ####
  
  
  # Model matrices:
  modMatrices <- list()
  
  # Obtain the expected homogeneous cor structure from means:
  expCorsVec <- model@sample@means[[1]]
  expCors <- matrix(1,nNode,nNode)
  expCors[lower.tri(expCors)] <- expCorsVec
  expCors[upper.tri(expCors)] <- t(expCors)[upper.tri(expCors)]
  
  
  # Fix sigma
  if (type == "cov"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'cov'. Use type = 'cor' or set corinput = FALSE")
    }
    
    modMatrices$sigma <- matrixsetup_sigma(sigma, 
                                           expcov=expCors,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = vars,
                                           equal = FALSE, 
                                           sampletable = sampleStats)    
  } else if (type == "chol"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'chol'.")
    }
    modMatrices$lowertri <- matrixsetup_lowertri(lowertri, 
                                                 expcov=expCors,
                                                 nNode = nNode, 
                                                 nGroup = nGroup, 
                                                 labels = vars,
                                                 equal = FALSE, 
                                                 sampletable = sampleStats)
  } else if (type == "ggm"){
    # Add omega matrix:
    modMatrices$omega <- matrixsetup_omega(omega, 
                                           expcov=expCors,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = vars,
                                           equal = FALSE, 
                                           sampletable = sampleStats)
    
    if (!corinput){
      # Add delta matrix (ingored if corinput == TRUE):
      modMatrices$delta <- matrixsetup_delta(delta, 
                                             expcov=expCors,
                                             nNode = nNode, 
                                             nGroup = nGroup, 
                                             labels = vars,
                                             equal = FALSE, 
                                             sampletable = sampleStats)       
    }
    
  } else if (type == "prec"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'prec'. Use type = 'ggm' or set corinput = FALSE")
    }
    
    # Add omega matrix:
    modMatrices$kappa <- matrixsetup_kappa(kappa, 
                                           expcov=expCors,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = vars,
                                           equal = FALSE, 
                                           sampletable = sampleStats)
  } else if (type == "cor"){
    # Add rho matrix:
    modMatrices$rho <- matrixsetup_rho(rho, 
                                       expcov=expCors,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = vars,
                                       equal = FALSE, 
                                       sampletable = sampleStats)
    
    if (!corinput){
      # Add SD matrix (ignored if corinput == TRUE):
      modMatrices$SD <- matrixsetup_SD(SD, 
                                       expcov=expCors,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = vars,
                                       equal = FALSE, 
                                       sampletable = sampleStats) 
    }      
  }
  
  
  # Compute expected random effects matrix:
  expRanEffects <- as.matrix(spectralshift(sampleStats@covs[[1]] - avgVmat))
  
  
  # Add random effets matrices:
  if (randomEffects == "cov"){
      modMatrices$sigma_randomEffects <- matrixsetup_sigma(sigma_randomEffects, 
                                           expcov=list(expRanEffects),
                                           nNode = nCor, 
                                           nGroup = nGroup, 
                                           labels = corvars,
                                           equal = FALSE, 
                                           sampletable = sampleStats, 
                                           name = "randomEffects")    
  } else if (randomEffects){
   
    modMatrices$lowertri_randomEffects <- matrixsetup_lowertri(lowertri_randomEffects, 
                                                               expcov=list(expRanEffects),
                                                               nNode = nCor, 
                                                               nGroup = nGroup, 
                                                               labels = corvars,
                                                 equal = FALSE, 
                                                 sampletable = sampleStats, 
                                                 name = "randomEffects")
  } else if (randomEffects == "ggm"){
    # Add omega matrix:
    modMatrices$omega_randomEffects <- matrixsetup_omega(omega_randomEffects, 
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor, 
                                                         nGroup = nGroup, 
                                                         labels = corvars,
                                           equal = FALSE, 
                                           sampletable = sampleStats, 
                                           name = "randomEffects")
    
   
      # Add delta matrix (ingored if corinput == TRUE):
      modMatrices$delta_randomEffects <- matrixsetup_delta(delta_randomEffects, 
                                                           expcov=list(expRanEffects),
                                                           nNode = nCor, 
                                                           nGroup = nGroup, 
                                                           labels = corvars,
                                             equal = FALSE, 
                                             sampletable = sampleStats, 
                                             name = "randomEffects")       

    
  } else if (randomEffects == "prec"){

    # Add omega matrix:
    modMatrices$kappa_randomEffects <- matrixsetup_kappa(kappa_randomEffects, 
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor, 
                                                         nGroup = nGroup, 
                                                         labels = corvars,
                                           equal = FALSE, 
                                           sampletable = sampleStats, 
                                           name = "randomEffects")
  } else if (randomEffects == "cor"){
    # Add rho matrix:
    modMatrices$rho_randomEffects <- matrixsetup_rho(rho_randomEffects, 
                                                     expcov=list(expRanEffects),
                                                     nNode = nCor, 
                                                     nGroup = nGroup, 
                                                     labels = corvars,
                                       equal = FALSE, 
                                       sampletable = sampleStats, 
                                       name = "randomEffects")
    

      # Add SD matrix (ignored if corinput == TRUE):
      modMatrices$SD_randomEffects <- matrixsetup_SD(SD_randomEffects, 
                                                     expcov=list(expRanEffects),
                                                     nNode = nCor, 
                                                     nGroup = nGroup, 
                                                     labels = corvars,
                                       equal = FALSE, 
                                       sampletable = sampleStats, 
                                       name = "randomEffects") 
  }
  
  
  ### UNTILLL HERE ###
  browser()
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  if (type == "cov" || type == "prec"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nCor), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nCor), # Elinimation matrix
      In = Diagonal(nCor) # Identity of dim n
    )
  } else if (type == "chol"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nCor), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nCor), # Elinimation matrix
      In = Diagonal(nCor), # Identity of dim n
      C = as(lavaan::lav_matrix_commutation(nCor,nCor),"sparseMatrix")
    )
  } else if (type == "ggm" || type == "cor"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nCor), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nCor), # Elinimation matrix
      Dstar = psychonetrics::duplicationMatrix(nCor,diag = FALSE), # Strict duplicaton matrix
      In = Diagonal(nCor), # Identity of dim n
      A = psychonetrics::diagonalizationMatrix(nCor)
    )
  }
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Form baseline model:
    if ("omega" %in% equal | "sigma" %in% equal | "kappa" %in% equal){
      equal <- c(equal,"lowertri")
    }
    if (!corinput){
      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "chol",
                                                  lowertri = "empty",
                                                  vars = corvars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  meanstructure=meanstructure,
                                                  corinput = corinput,
                                                  ordered = ordered,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats)
    } else {
      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "cor",
                                                  lowertri = "empty",
                                                  vars = corvars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  meanstructure=meanstructure,
                                                  corinput = corinput,
                                                  ordered = ordered,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats)
    }
    
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    if (!corinput){
      model@baseline_saturated$saturated <- varcov(data,
                                                   type = "chol", 
                                                   lowertri = "full", 
                                                   vars = vars,
                                                   groups = groups,
                                                   covs = covs,
                                                   means = means,
                                                   nobs = nobs,
                                                   missing = missing,
                                                   equal = equal,
                                                   estimator = estimator,
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats)
    } else {
      model@baseline_saturated$saturated <- varcov(data,
                                                   type = "cor", 
                                                   lowertri = "full", 
                                                   vars = vars,
                                                   groups = groups,
                                                   covs = covs,
                                                   means = means,
                                                   nobs = nobs,
                                                   missing = missing,
                                                   equal = equal,
                                                   estimator = estimator,
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats)
    }
    
    
    # if not FIML, Treat as computed:
    if (estimator != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)      
    }
  }
  
  # Return model:
  return(model)
}