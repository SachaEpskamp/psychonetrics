# ggm model creator:
ggm <- function(
  data, # Dataset
  adjacency = "empty", # (only lower tri is used) "empty", "full" or adjacency structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml",
  equal = "none", # Can also be: c("network","means")
  mu
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
  model <- generate_psychonetrics(model = "ggm",sample = sampleStats,computed = FALSE)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Fix the adjacency matrix:
  adjacency <- fixAdj(adjacency,nGroup,nNode,"network" %in% equal)
  
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
    list(adjacency,
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
  model@fitfunctions <- list(
    fitfunction = fit_ggm,
    gradient = gradient_ggm,
    hessian = hessian_ggm,
    extramatrices = list(
        D = as(matrixcalc::duplication.matrix(nNode),"sparseMatrix"),
        M = Mmatrix(model@parameters)
      )
  )
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  

  # Return model:
  return(model)
}