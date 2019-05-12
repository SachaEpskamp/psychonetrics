# Latent network model creator
varcov <- function(
  data, # Dataset
  type = c("cov","chol","prec","ggm"), # Maybe add cor at some point, but not now
  sigma = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa = "full", # Precision
  # rho = "full", # Correlations
  omega = "full", # Partial correlations
  lowertri = "full", # Cholesky
  delta = "full", # Used for both ggm and pcor
  mu,
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "ML",
  optimizer = "default",
  storedata = FALSE,
  WLS.V,
  sampleStats # Leave to missing
){
  rawts = FALSE
  if (rawts){
    warning("'rawts' is only included for testing purposes. Please do not use!")
  }
  
  # Type:
  type <- match.arg(type)

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
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing = ifelse(estimator == "FIML","pairwise",missing),
                               rawts = rawts,
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               weightsmatrix = WLS.V)
  }


  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "varcov",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(y = type),
                                  submodel = type)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode+1) / 2 * nGroup + # Covariances per group
    nNode * nGroup # Means per group
  
  # Model matrices:
  modMatrices <- list()
  
  # Fix mu
  modMatrices$mu <- matrixsetup_mu(mu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "mu" %in% equal,
                       expmeans = model@sample@means, sampletable = sampleStats)
  
  # fixMu(mu,nGroup,nNode,"mu" %in% equal)
  
  # Fix sigma
  if (type == "cov"){
    modMatrices$sigma <- matrixsetup_sigma(sigma, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "sigma" %in% equal, sampletable = sampleStats)    
  } else if (type == "chol"){
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
    
    # Add delta matrix:
    modMatrices$delta <- matrixsetup_delta(delta, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "delta" %in% equal, sampletable = sampleStats) 
  } else if (type == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa <- matrixsetup_kappa(kappa, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "kappa" %in% equal, sampletable = sampleStats)
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
  } else if (type == "ggm"){
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
                                             baseline_saturated = FALSE,sampleStats=sampleStats)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
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
           baseline_saturated = FALSE,sampleStats=sampleStats)
    
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