# Latent network model creator
varcov <- function(
  data, # Dataset
  type = c("cov","chol","cor","prec","pcor"),
  sigma = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa = "full", # Precision
  rho = "full", # Correlations
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
  optimizer = "default"
){
  rawts = FALSE
  if (rawts){
    warning("'rawts' is only included for testing purposes. Please do not use!")
  }
  
  # Type:
  type <- match.arg(type)

  # Obtain sample stats:
  sampleStats <- samplestats(data = data, 
                             vars = vars, 
                             groups = groups,
                             covs = covs, 
                             means = means, 
                             nobs = nobs, 
                             missing = ifelse(estimator == "FIML","pairwise",missing),
                             rawts = rawts,
                             fimldata = estimator == "FIML")

  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "varcov",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts)
  
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
  modMatrices$sigma <- matrixsetup_sigma(sigma, 
                                  expcov=model@sample@covs,
                nNode = nNode, 
                nGroup = nGroup, 
                labels = sampleStats@variables$label,
                equal = "sigma" %in% equal, sampletable = sampleStats)
  
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  model@extramatrices <- list(
    D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    In = Diagonal(nNode) # Identity of dim n
  )
    
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
   
    # Form baseline model:
    model@baseline_saturated$baseline <- cholesky(data, 
                                             lowertri = "empty",
                                             vars = vars,
                                             groups = groups,
                                             covs = covs,
                                             means = means,
                                             nobs = nobs,
                                             missing = missing,
                                             equal = equal,
                                             estimator = estimator,
                                             baseline_saturated = FALSE)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- cholesky(data, 
           lowertri = "full", 
           vars = vars,
           groups = groups,
           covs = covs,
           means = means,
           nobs = nobs,
           missing = missing,
           equal = equal,
           estimator = estimator,
           baseline_saturated = FALSE)
    
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