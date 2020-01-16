# Latent network model creator
varcov <- function(
  data, # Dataset
  type = c("cov","chol","prec","ggm","cor"),
  sigma = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa = "full", # Precision
  # rho = "full", # Correlations
  omega = "full", # Partial correlations
  lowertri = "full", # Cholesky
  delta = "full", # Used for both ggm and pcor
  rho = "full", # Used for cor
  SD = "full", # Used for cor
  mu,
  tau,
  vars, # character indicating the variables Extracted if missing from data - group variable
  ordered = character(0), # character indicating the variables that are ordinal
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer = "default",
  storedata = FALSE,
  WLS.V,
  sampleStats, # Leave to missing
  meanstructure, # Defaults to TRUE if data is used or means is used, FALSE otherwie
  corinput, # Defaults to TRUE if the input is detected to consist of correlation matrix/matrices, FALSE otherwise
  verbose = TRUE,
  covtype = c("choose","ML","UB")
){
  rawts = FALSE
  if (rawts){
    warning("'rawts' is only included for testing purposes. Please do not use!")
  }
  
  # Reset ordered if needed:
  if (identical(ordered, FALSE)){
    ordered <- character(0)
  }
  
  if (estimator == "default"){
    if (length(ordered) > 0){
      estimator <- "WLS"
    } else {
      estimator <- "ML"
    }
  }
  
  # Check WLS for ordinal:
  if (length(ordered) > 0 & !estimator %in% c("WLS","DWLS","ULS")){
    stop("Ordinal data is only supported for WLS, DWLS and ULS estimators.")
  }
  
  # Check corinput FIXME: may be mixture
  if (length(ordered) > 0 & missing(corinput)){
    corinput <- TRUE
  }
  if (length(ordered) > 0 && (!missing(corinput) && !corinput)){
    stop("corinput must be TRUE for ordinal data.")
  }

  
  # Disable meanstructure (FIXME: Allow for both ordinal and continuous):
  if (length(ordered) > 0){
    meanstructure <- FALSE
  }
  
  # Type:
  type <- match.arg(type)
  
  # Set meanstructure:
  if (missing(meanstructure)){
    meanstructure <- (!missing(data) || !missing(means))
  }
  
  # Check FIML:
  if (!meanstructure && estimator == "FIML"){
    stop("meanstructure = FALSE is not yet supported for 'FIML' estimator")
  }
  
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    # WLS weights:
    if (missing(WLS.V)){
      WLS.V <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                      switch(estimator,
                             "WLS" = "full",
                             "ULS" = "identity",
                             "DWLS" = "diag"
                      ))
    }
    
    sampleStats <- samplestats(data = data, 
                               vars = vars, 
                               ordered=ordered,
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing = ifelse(estimator == "FIML","pairwise",missing),
                               rawts = rawts,
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               weightsmatrix = WLS.V,
                               meanstructure = meanstructure,
                               corinput = corinput,
                               covtype=covtype,
                               verbose=verbose)
  }
 
  # Overwrite corinput:
  corinput <- sampleStats@corinput
  
  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "varcov",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(y = type),
                                  submodel = type, meanstructure = meanstructure)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  

  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  if (length(ordered) > 0){
    nThresh <- sum(sapply(model@sample@thresholds,function(x)sum(sapply(x,length))))    
  } else {
    nThresh <- 0
  }

  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode-1) / 2 * nGroup + # Covariances per group
    (!corinput) * nNode * nGroup + # Variances (ignored if correlation matrix is input)
    meanstructure * nMeans + # Means per group
    nThresh
  
  # Model matrices:
  modMatrices <- list()
  
  # # Fix mu
  # modMatrices$mu <- matrixsetup_mu(mu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "mu" %in% equal,
  #                                  expmeans = model@sample@means, sampletable = sampleStats, meanstructure = meanstructure)
  
  # Alternative, without mean structure just ignore mu:
  if (meanstructure){
    # Fix mu
    modMatrices$mu <- matrixsetup_mu(mu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "mu" %in% equal,
                       expmeans = model@sample@means, sampletable = sampleStats, meanstructure = meanstructure)
  }
  
  # Thresholds:
  if (length(ordered) > 0){
    # Ideal setup for tau is a matrix, with NA for the missing elements.
    modMatrices$tau <- matrixsetup_tau(tau, nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,
                                       equal = "tau" %in% equal, sampleThresholds = model@sample@thresholds, sampletable = sampleStats)
    
  }
  
  # fixMu(mu,nGroup,nNode,"mu" %in% equal)
  
  # Fix sigma
  if (type == "cov"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'cov'. Use type = 'cor' or set corinput = FALSE")
    }
    
    modMatrices$sigma <- matrixsetup_sigma(sigma, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "sigma" %in% equal, sampletable = sampleStats)    
  } else if (type == "chol"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'chol'.")
    }
    modMatrices$lowertri <- matrixsetup_lowertri(lowertri, 
                                                 expcov=model@sample@covs,
                                                 nNode = nNode, 
                                                 nGroup = nGroup, 
                                                 labels = sampleStats@variables$label,
                                                 equal = "lowertri" %in% equal, sampletable = sampleStats)
  } else if (type == "ggm"){
    # Add omega matrix:
    modMatrices$omega <- matrixsetup_omega(omega, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "omega" %in% equal, sampletable = sampleStats)
    
    if (!corinput){
      # Add delta matrix (ingored if corinput == TRUE):
      modMatrices$delta <- matrixsetup_delta(delta, 
                                             expcov=model@sample@covs,
                                             nNode = nNode, 
                                             nGroup = nGroup, 
                                             labels = sampleStats@variables$label,
                                             equal = "delta" %in% equal, sampletable = sampleStats)       
    }
    
  } else if (type == "prec"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'prec'. Use type = 'ggm' or set corinput = FALSE")
    }
    
    # Add omega matrix:
    modMatrices$kappa <- matrixsetup_kappa(kappa, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "kappa" %in% equal, sampletable = sampleStats)
  } else if (type == "cor"){
    # Add rho matrix:
    modMatrices$rho <- matrixsetup_rho(rho, 
                                       expcov=model@sample@covs,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = sampleStats@variables$label,
                                       equal = "rho" %in% equal, sampletable = sampleStats)
    
    if (!corinput){
      # Add SD matrix (ignored if corinput == TRUE):
      modMatrices$SD <- matrixsetup_SD(SD, 
                                       expcov=model@sample@covs,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = sampleStats@variables$label,
                                       equal = "SD" %in% equal, sampletable = sampleStats) 
    }      
  }
  
  
  
  
  
 
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  if (type == "cov" || type == "prec"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
      In = Diagonal(nNode) # Identity of dim n
    )
  } else if (type == "chol"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
      In = Diagonal(nNode), # Identity of dim n
      C = as(lavaan::lav_matrix_commutation(nNode,nNode),"sparseMatrix")
    )
  } else if (type == "ggm" || type == "cor"){
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
      Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
      In = Diagonal(nNode), # Identity of dim n
      A = psychonetrics::diagonalizationMatrix(nNode)
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
      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "cor",
                                                  lowertri = "empty",
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