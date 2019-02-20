# precision model creator:
precision <- function(
  data, # Dataset
  kappa = "empty", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml",
  equal = "none", # Can also be: c("network","means")
  mu,
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  fitfunctions # Leave empty
){
  # Obtain sample stats:
  sampleStats <- samplestats(data = data, 
                             vars = vars, 
                             groups = groups,
                             covs = covs, 
                             means = means, 
                             nobs = nobs, 
                             missing = missing)
  
  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "precision",sample = sampleStats,computed = FALSE, equal = equal)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode+1) / 2 * nGroup + # Covariances per group
    nNode * nGroup # Means per group
  
  # Fix the kappa matrix:
  kappa <- fixAdj(kappa,nGroup,nNode,"network" %in% equal)
  
  # Check mu?
  mu <- fixMu(mu,nGroup,nNode,"means" %in% equal)

  # Generate the full parameter table:
  pars <- generateAllParameterTables(
    # Mu:
    list(mu,
         mat =  "mu",
         op =  "~1",
         symmetrical= FALSE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = sampleStats@variables$label),
    
    # Kappa:
    list(kappa,
         mat =  "kappa",
         op =  "--",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = sampleStats@variables$label,
         sparse = TRUE,
         posdef = TRUE
      )
  )
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  # Form the fitfunctions list:
  if (missing(fitfunctions)){
    model@fitfunctions <- list(
      fitfunction = fit_precision,
      gradient = gradient_precision,
      hessian = hessian_precision,
      loglik=loglik_precision
      # extramatrices = list(
      #   D = as(matrixcalc::duplication.matrix(nNode),"sparseMatrix"),
      #   M = Mmatrix(model@parameters)
      # )
    )    
  } else {
    model@fitfunctions <- fitfunctions
  }

    
    # Form the model matrices
    model@modelmatrices <- formModelMatrices(model)

  
  ### Baseline model ###
  if (baseline_saturated){
    model@baseline_saturated$baseline <- precision(data = data, 
                                             kappa = "empty",
                                             vars = vars,
                                             groups = groups,
                                             covs = covs,
                                             means = means,
                                             nobs = nobs,
                                             missing = missing,
                                             equal = equal,
                                             baseline_saturated = FALSE,
                                             fitfunctions = model@fitfunctions) 
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    # Run:
    # model@baseline_saturated$baseline <- runmodel(model@baseline_saturated$baseline, addfit = FALSE, addMIs = FALSE)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- precision(data = data, 
                                              kappa = "full",
                                              vars = vars,
                                              groups = groups,
                                              covs = covs,
                                              means = means,
                                              nobs = nobs,
                                              missing = missing,
                                              equal = "none",
                                              baseline_saturated = FALSE,
                                              fitfunctions = model@fitfunctions)
    
    # Add model:
    # model@baseline_saturated$saturated@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$saturated@parameters)
    
    # Run:
    # model@baseline_saturated$saturated <- runmodel(model@baseline_saturated$saturated, addfit = FALSE, addMIs = FALSE)
  }

  
  # Return model:
  return(model)
}