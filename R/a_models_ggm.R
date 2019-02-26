# Latent network model creator
ggm <- function(
  data, # Dataset
  omega = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta = "full", # If missing, just full for both groups or equal
  mu,
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "ML",
  optimizer = "ucminf"
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
  model <- generate_psychonetrics(model = "ggm",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode+1) / 2 * nGroup + # Covariances per group
    nNode * nGroup # Means per group
  
  # Fix mu
  mu <- fixMu(mu,nGroup,nNode,"mu" %in% equal)
  
  # Fix omega
  omega <- fixAdj(omega,nGroup,nNode,"omega" %in% equal,diag0=TRUE)
  
  # Fix delta:
  delta <- fixAdj(delta,nGroup,nNode,"delta" %in% equal,diagonal=TRUE)
  
  
  # Generate starting values:
  omegaStart <- omega
  deltaStart <- delta
  muStart <- mu
  for (g in 1:nGroup){
    # Means with sample means:
    muStart[,g] <- sampleStats@means[[g]]
    
    # For the network, compute a rough glasso network:
    zeroes <- which(omegaStart[,,g]==0 & t(omegaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      wi <- corpcor::pseudoinverse(sampleStats@covs[[g]])
    } else {
      glas <- glasso(as.matrix(sampleStats@covs[[g]]),
                     rho = 1e-10, zero = zeroes)
      wi <- glas$wi
    }

    # Network starting values:
    omegaStart[,,g] <- as.matrix(qgraph::wi2net(wi))
    diag(omegaStart[,,g] ) <- 0
    
    # Delta:
    deltaStart[,,g] <- diag(1/sqrt(diag(wi)))
  }

  # Generate the full parameter table:
  pars <- generateAllParameterTables(
    # Mu:
    list(mu,
         mat =  "mu",
         op =  "~1",
         symmetrical= FALSE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = "1",
         start = muStart),
    
 
    # Omega:
    list(omega,
         mat =  "omega",
         op =  "--",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = sampleStats@variables$label,
         sparse = TRUE,
         posdef = TRUE,
         diag0=TRUE,
         lower = -1,
         upper = 1,
         start = omegaStart
    ),
    
    # Delta_eta:
    list(delta,
         mat =  "delta",
         op =  "~/~",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = sampleStats@variables$label,
         sparse = TRUE,
         posdef = TRUE,
         diagonal = TRUE,
         lower = 0,
         start = deltaStart
    )
    
  )
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  model@extramatrices <- list(
    D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    In = Diagonal(nNode), # Identity of dim n
    A = psychonetrics::diagonalizationMatrix(nNode)
  )
    
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
   
    # Form baseline model:
    model@baseline_saturated$baseline <- ggm(data, 
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
    model@baseline_saturated$saturated <- ggm(data, 
           omega = "full", 
           vars = vars,
           groups = groups,
           covs = covs,
           means = means,
           nobs = nobs,
           missing = missing,
           equal = equal,
           baseline_saturated = FALSE)
    
    # Treat as computed:
    model@baseline_saturated$saturated@computed <- TRUE
    
    # FIXME: TODO
    model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)

  }
  
  # Return model:
  return(model)
}