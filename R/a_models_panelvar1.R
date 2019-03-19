# Latent network model creator
panelvar1 <- function(
  data, # Dataset
  vars, # Must be a matrix!
  
  # Type:
  contemporaneous = c("cov","chol","prec","ggm"), 
  between = c("cov","chol","prec","ggm"), 
  
  # Temporal effects:
  beta = "full",
  
  # Contemporaneous effects:
  omega_zeta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_zeta = "full", # If missing, just full for both groups or equal
  kappa_zeta = "full",
  sigma_zeta = "full",
  lowertri_zeta = "full",
  
  # Between-case effects:
  omega_mu = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_mu = "full", # If missing, just full for both groups or equal
  kappa_mu = "full",
  sigma_mu = "full",
  lowertri_mu = "full",
  
  # The rest:
  mu_y, # Grand mean
  beepvar,
  dayvar,
  idvar,
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  # fitfunctions, # Leave empty
  estimator = "ML",
  optimizer = "default"
){
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  
  # Extract var names:
  if (missing(vars)){
    stop("'vars' argument may not be missing")
  }
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix, with rows indicating variables and columns indicating measurements.")
  }
  
  # List all variables to use, in order:
  allVars <- na.omit(as.vector(vars))
  
  # Obtain sample stats:
  sampleStats <- samplestats(data = data, 
                             vars = allVars, 
                             groups = groups,
                             covs = covs, 
                             means = means, 
                             nobs = nobs, 
                             missing  = ifelse(estimator == "FIML","pairwise",missing),
                             fimldata = estimator == "FIML")
  
  # Design matrix:
  design <- as(1*(!is.na(vars)),"sparseMatrix")
  
  # time per var:
  timePerVar <- as.vector(design * row(design))
  timePerVar <- timePerVar[timePerVar!=0]
  
  # Number of nodes:
  nNode <- nrow(vars)
  
  # Number of measurements:
  nTime <- ncol(vars)
  
  # row names:
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("V",seq_len(nrow(vars)))
  }
  varnames <- rownames(vars)
  
  # col names:
  if (is.null(colnames(vars))){
    colnames(vars) <- paste0("T",seq_len(ncol(vars)))
  }
  timenames <- colnames(vars)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "panelvar1", types = list(contemporaneous = contemporaneous, between = between),
                                  sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # FIXME: Keep this the same for now for rawts = TRUE
  nVar <- length(allVars)
  
  # Add number of observations:
  model@sample@nobs <-  
    nVar * (nVar+1) / 2 * nGroup + # Covariances per group
    nVar * nGroup # Means per group
  
  # Model matrices:
  modMatrices <- list()
  
  # Expected means:
  expMeans <- lapply(model@sample@means, function(x)tapply(x,timePerVar,mean,na.rm=TRUE))
  
  # Fix mu_y
  modMatrices$mu_y <- matrixsetup_mu(mu_y,nNode = nNode, nGroup = nGroup, labels = varnames,equal = "mu_y" %in% equal,
                                   expmeans = expMeans, sampletable = sampleStats, name = "mu_y")
  
  
  # FIXME: simple start values for now...
  
  # T=1 cov structure:
  firstVars <- apply(vars,1,function(x)na.omit(x)[1])
  firstCovs <- lapply(sampleStats@covs,function(x)spectralshift(x[firstVars,firstVars])/3)
  
  # Setup Beta:
  modMatrices$beta <- matrixsetup_beta(beta, 
                                       name = "beta",
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = varnames,
                                       equal = "beta" %in% equal, sampletable = sampleStats)
  
  
  # Between-case effects:
  if (contemporaneous == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta, name = "sigma_zeta",
                                                expcov=firstCovs,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = varnames,
                                                equal = "sigma_zeta" %in% equal, sampletable = sampleStats)    
  } else if (contemporaneous == "chol"){
    modMatrices$lowertri_zeta <- matrixsetup_lowertri(lowertri_zeta,  name = "lowertri_zeta",
                                                      expcov=firstCovs,
                                                      nNode = nNode, 
                                                      nGroup = nGroup, 
                                                      labels = varnames,
                                                      equal = "lowertri_zeta" %in% equal, sampletable = sampleStats)
  } else if (contemporaneous == "ggm"){
    # Add omega matrix:
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta, name = "omega_zeta",
                                                expcov=firstCovs,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = varnames,
                                                equal = "omega_zeta" %in% equal, sampletable = sampleStats)
    
    # Add delta matrix:
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta, name = "delta_zeta",
                                                expcov=firstCovs,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = varnames,
                                                equal = "delta_zeta" %in% equal, sampletable = sampleStats) 
  } else if (contemporaneous == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_zeta <- matrixsetup_kappa(kappa_zeta, name = "kappa_zeta",
                                                expcov=firstCovs,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = varnames,
                                                equal = "kappa_zeta" %in% equal, sampletable = sampleStats)
  }
  
  # Between-case effects:
  if (between == "cov"){
    modMatrices$sigma_mu <- matrixsetup_sigma(sigma_mu, name = "sigma_mu",
                                              expcov=firstCovs,
                                              nNode = nNode, 
                                              nGroup = nGroup, 
                                              labels = varnames,
                                              equal = "sigma_mu" %in% equal, sampletable = sampleStats)    
  } else if (between == "chol"){
    modMatrices$lowertri_mu <- matrixsetup_lowertri(lowertri_mu,  name = "lowertri_mu",
                                                    expcov=firstCovs,
                                                    nNode = nNode, 
                                                    nGroup = nGroup, 
                                                    labels = varnames,
                                                    equal = "lowertri_mu" %in% equal, sampletable = sampleStats)
  } else if (between == "ggm"){
    # Add omega matrix:
    modMatrices$omega_mu <- matrixsetup_omega(omega_mu, name = "omega_mu",
                                              expcov=firstCovs,
                                              nNode = nNode, 
                                              nGroup = nGroup, 
                                              labels = varnames,
                                              equal = "omega_mu" %in% equal, sampletable = sampleStats)
    
    # Add delta matrix:
    modMatrices$delta_mu <- matrixsetup_delta(delta_mu, name = "delta_mu",
                                              expcov=firstCovs,
                                              nNode = nNode, 
                                              nGroup = nGroup, 
                                              labels = varnames,
                                              equal = "delta_mu" %in% equal, sampletable = sampleStats) 
  } else if (between == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_mu <- matrixsetup_kappa(kappa_mu, name = "kappa_mu",
                                              expcov=firstCovs,
                                              nNode = nNode, 
                                              nGroup = nGroup, 
                                              labels = varnames,
                                              equal = "kappa_mu" %in% equal, sampletable = sampleStats)
  }

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  model@extramatrices <- list(
    D =  psychonetrics::duplicationMatrix(nVar), # Toeplitz matrix D 
    D2 = psychonetrics::duplicationMatrix(nNode), # Toeplitz matrix D 
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    In = Diagonal(nNode), # Identity of dim n
    A = psychonetrics::diagonalizationMatrix(nNode),
    C = as(lavaan::lav_matrix_commutation(nNode,nNode),"sparseMatrix"),
    design = design
    # P=P # Permutation matrix
  )
  
  
  # Come up with P...
  # Dummy matrix to contain indices:
  # Dummy matrices with indices:
  muDummy <- matrix(rep(1:nNode,nTime))
  sigDummy <- matrix(0,nNode,nNode)
  sigDummy[lower.tri(sigDummy,diag=TRUE)] <- max(muDummy) + seq_len(nNode*(nNode+1)/2)
  sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]
  
  U <- list(sigDummy)
  # Now make all lag-k blocks...
  # Form blocks:
  for (i in 1:(nTime-1)){
    U[[length(U) + 1]] <- matrix(max(unlist(U)) + seq_len(nNode^2), nNode, nNode)
  }
  
  allSigmas <- blockToeplitz(U)
  
  # Now subset with only observed:
  subMu <- muDummy[as.vector(design==1),,drop=FALSE]
  subSigmas <- allSigmas[as.vector(design==1),as.vector(design==1)]
  
  inds <- c(as.vector(subMu),subSigmas[lower.tri(subSigmas,diag=TRUE)])
  
  # P matrix:
  # P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=inds))
  distVec <-  c(as.vector(subMu),subSigmas[lower.tri(subSigmas,diag=TRUE)])
  nTotal <- length(distVec)
  distVecrawts <- seq_along(distVec)[distVec!=0]
  distVec <- distVec[distVec!=0]
  # Now I can make the matrix:
  
  model@extramatrices$P <- sparseMatrix(
    i = distVecrawts, j = distVec, dims = c(nTotal, max(subSigmas))
  )
  
  # model@extramatrices$P <- sparseMatrix(j=seq_along(inds),i=order(inds))
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Form baseline model:
    model@baseline_saturated$baseline <- varcov(data,
                                                type = "chol",
                                                lowertri = "empty",
                                                vars = allVars,
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
    model@baseline_saturated$saturated <- varcov(data,
                                                 type = "chol", 
                                                 lowertri = "full", 
                                                 vars = allVars,
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
