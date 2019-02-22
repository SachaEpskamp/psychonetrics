# Latent network model creator
lnm <- function(
  data, # Dataset
  lambda, # only required non-missing matrix
  omega_eta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_eta = "full", # If missing, just full for both groups or equal
  sigma_epsilon = "diag", # Residuals, will default to a diagonal matrix
  tau,
  mu_eta,
  vars, # character indicating the variables Extracted if missing from data - group variable
  latents, # Name of latent varianles
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  fitfunctions, # Leave empty
  identify = TRUE
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
  model <- generate_psychonetrics(model = "lnm",sample = sampleStats,computed = FALSE, equal = equal)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode+1) / 2 * nGroup + # Covariances per group
    nNode * nGroup # Means per group
  
  # Stop if lambda is missing or a character:
  if (missing(lambda)){
    stop("'lambda' may not be missing")
  }
  if (is.character(lambda)){
    stop("'lambda' may not be a string")
  }
  
  # Fix tau
  tau <- fixMu(tau,nGroup,nNode,"tau" %in% equal)
  
  # Fix the lambda matrix:
  lambda <- fixMatrix(lambda,nGroup,equal = "lambda" %in% equal)
  
  # Number of latents:
  nLatent <- dim(lambda)[2]
  
  # Fix mu_eta
  mu_eta <- fixMu(mu_eta,nGroup,nLatent,"tau" %in% equal)
  
  # Fix the omega_eta matrix:
  omega_eta <- fixAdj(omega_eta,nGroup,nLatent,"omega_eta" %in% equal,diag0=TRUE)
  
  # Fix delta:
  delta_eta <- fixAdj(delta_eta,nGroup,nLatent,"delta_eta" %in% equal,diagonal=TRUE)
  
  # Fix sigma_epsilon
  if (is.character(sigma_epsilon) && sigma_epsilon == "empty"){
    diag0 <- TRUE
  } else {
    diag0 <- FALSE
  }
  if (is.character(sigma_epsilon) && sigma_epsilon == "diag"){
    sigma_epsilon <- "empty"
  }
  sigma_epsilon <- fixAdj(sigma_epsilon,nGroup,nNode,"sigma_epsilon" %in% equal,diag0=diag0)
  
  # If latents is not provided, make it:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }
  
  
  # Generate the full parameter table:
  pars <- generateAllParameterTables(
    # Tau:
    list(tau,
         mat =  "tau",
         op =  "~1",
         symmetrical= FALSE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = "1"),
    
    # mu_eta:
    list(mu_eta,
         mat =  "mu_eta",
         op =  "~1",
         symmetrical= FALSE, 
         sampletable=sampleStats,
         rownames = latents,
         colnames = "1"),
    
    # Lambda:
    list(lambda,
         mat =  "lambda",
         op =  "~=",
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = latents,
         sparse = TRUE
    ),
    
    # Omega:
    list(omega_eta,
         mat =  "omega_eta",
         op =  "--",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = latents,
         colnames = latents,
         sparse = TRUE,
         posdef = TRUE,
         diag0=TRUE,
         lower = -1,
         upper = 1
    ),
    
    # Delta_eta:
    list(delta_eta,
         mat =  "delta_eta",
         op =  "~/~",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = latents,
         colnames = latents,
         sparse = TRUE,
         posdef = TRUE,
         diagonal = TRUE,
         lower = 1e-10
    ),
    
    # Sigma_epsilon:
    list(sigma_epsilon,
         mat =  "sigma_epsilon",
         op =  "~~",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = sampleStats@variables$label,
         sparse = TRUE,
         posdef = TRUE,
         lower = 1e-10
    )
    
  )
  
  # Set lambda start:
  pars$partable$est[pars$partable$matrix=="lambda" & !pars$partable$fixed] <- 0.5

  # Set tau startL
  for (g in 1:nGroup){
    pars$partable$est[pars$partable$matrix=="tau" & pars$partable$group_id == g] <- sampleStats@means[[g]]
  }
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  # Form the fitfunctions list:
  if (missing(fitfunctions)){
    model@fitfunctions <- list(
      fitfunction = fit_lnm,
      gradient = gradient_lnm,
      # hessian = hessian_ggm,
      loglik=loglik_ggm,
      information = Fisher_lnm,
      extramatrices = list(
        D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
        L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
        Dstar = psychonetrics::duplicationMatrix(nLatent,diag = FALSE), # Strict duplicaton matrix
        ITheta = Diagonal(nNode*(nNode+1)/2),
        # A = diagonalizationMatrix(nNode), # Diagonalization matrix
        # An2 = diagonalizationMatrix(nNode^2), # Larger diagonalization matrix
        In = Diagonal(nNode), # Identity of dim n
        Inlatent = Diagonal(nLatent),
        C = as(lavaan::lav_matrix_commutation(nNode, nLatent),"sparseMatrix"),
        A = psychonetrics::diagonalizationMatrix(nLatent)
        # In2 = Diagonal(nNode^2), # Identity of dim n^2
        # In3 = Diagonal(nNode^3), # Identity of dim n^3
        # E = basisMatrix(nNode) # Basis matrix
      )
    )    
  } else {
    model@fitfunctions <- fitfunctions
  }
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    model@baseline_saturated$baseline <- ggm(data = data, 
                                             omega = "empty",
                                             vars = vars,
                                             groups = groups,
                                             covs = covs,
                                             means = means,
                                             nobs = nobs,
                                             missing = missing,
                                             equal = equal,
                                             baseline_saturated = FALSE) 
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- ggm(data = data, 
                                              omega = "full",
                                              vars = vars,
                                              groups = groups,
                                              covs = covs,
                                              means = means,
                                              nobs = nobs,
                                              missing = missing,
                                              equal = "none",
                                              baseline_saturated = FALSE)
    
    # Add model:
    # model@baseline_saturated$saturated@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$saturated@parameters)
    
    # Run:
    # model@baseline_saturated$saturated <- runmodel(model@baseline_saturated$saturated, addfit = FALSE, addMIs = FALSE)
  }
  
  # Identify model:
  if (identify){
    model <- identify(model)
  }
  
  # Return model:
  return(model)
}