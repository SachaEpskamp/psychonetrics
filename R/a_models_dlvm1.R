# Latent network model creator
dlvm1 <- function(
  data, # Dataset
  vars, # Must be a matrix!
  
  # Factor loadings:
  lambda, # May not be missing
  # lambda_within, # May not be missing
  # lambda_between, # May not be missing
  
  # Type:
  within_latent = c("cov","chol","prec","ggm"), 
  within_residual = c("cov","chol","prec","ggm"), 
  between_latent = c("cov","chol","prec","ggm"), 
  between_residual = c("cov","chol","prec","ggm"), 
  
  # Temporal effects:
  beta = "full",
  
  # Contemporaneous latent effects within:
  omega_zeta_within = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_zeta_within = "full", # If missing, just full for both groups or equal
  kappa_zeta_within = "full",
  sigma_zeta_within = "full",
  lowertri_zeta_within = "full",
  
  # Residual latent effects within:
  omega_epsilon_within = "empty", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_epsilon_within = "empty", # If missing, just full for both groups or equal
  kappa_epsilon_within = "empty",
  sigma_epsilon_within = "empty",
  lowertri_epsilon_within = "empty",
  
  # Contemporaneous latent effects between:
  omega_zeta_between = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_zeta_between = "full", # If missing, just full for both groups or equal
  kappa_zeta_between = "full",
  sigma_zeta_between = "full",
  lowertri_zeta_between = "full",
  
  # Residual latent effects between:
  omega_epsilon_between = "empty", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_epsilon_between = "empty", # If missing, just full for both groups or equal
  kappa_epsilon_between = "empty",
  sigma_epsilon_between = "empty",
  lowertri_epsilon_between = "empty",
  
  # Mean structure:
  nu,
  mu_eta,
  
  # Identification:
  identify = TRUE,
  identification = c("loadings","variance"),
  
  # Latents:
  # latents,
  # latents,
  latents,
  
  # The rest:
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  covtype = c("choose","ML","UB"),
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  # fitfunctions, # Leave empty
  estimator = "ML",
  optimizer,
  storedata = FALSE,
  verbose = FALSE,
  sampleStats
){
  covtype <- match.arg(covtype)
  # Check for missing:
  # if (missing(lambda_within)){
  #   stop("'lambda_within' may not be missing")
  # }
  # if (is.character(lambda_within)){
  #   stop("'lambda_within' may not be a string")
  # }
  # if (missing(lambda_between)){
  #   stop("'lambda_between' may not be missing")
  # }
  # if (is.character(lambda_between)){
  #   stop("'lambda_between' may not be a string")
  # }
  # 
  # if (missing(lambda)){
  #   stop("'lambda' may not be missing")
  # }
  # if (is.character(lambda)){
  #   stop("'lambda' may not be a string")
  # }
  
  
  # Match args:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  within_residual <- match.arg(within_residual)
  between_residual <- match.arg(between_residual)
  identification <- match.arg(identification)
  
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
  if (missing(sampleStats)){
    sampleStats <- samplestats(data = data, 
                               vars = allVars, 
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing  = ifelse(estimator == "FIML","pairwise",missing),
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               covtype=covtype,
                               verbose=verbose,
                               weightsmatrix = ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                                                      switch(estimator,
                                                        "WLS" = "full",
                                                        "ULS" = "identity",
                                                        "DWLS" = "diag"
                                                      )))
  }
  
  
  # Design matrix:
  # design <- as(1*(!is.na(vars)),"dgCMatrix")
  design <- as.matrix(1*(!is.na(vars)))
  
  # time per var:
  timePerVar <- as.vector(design * row(design))
  timePerVar <- timePerVar[timePerVar!=0]
  
  # Number of variables:
  nVar <- nrow(vars)
  
  # Check lambda:
  if (missing(lambda)){
    if (verbose){
      message("'lambda' is missing, creating observed data only model.")
    }
    
    lambda <- diag(nVar)
    O <- matrix(0, nVar, nVar)
    omega_epsilon_within <- O
    delta_epsilon_within <- O
    kappa_epsilon_within <- O
    sigma_epsilon_within <- O
    lowertri_epsilon_within <- O
    
    omega_epsilon_between <- O
    delta_epsilon_between <- O
    kappa_epsilon_between <- O
    sigma_epsilon_between <- O
    lowertri_epsilon_between <- O
  }
  
  # Number of measurements:
  nTime <- ncol(vars)
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  # Number of latents at within level:
  # nLat <- ncol(lambda_within)
  
  # Number of latents at between level:
  # nLat <- ncol(lambda_between)
  
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
  
  # Latents:
  # if (missing(latents)){
  #   latents <- paste0("Eta_within_",seq_len(nLat))
  # }
  # if (length(latents) != nLat){
  #   stop("Length of 'latents' is not equal to number of latent variables in model.")
  # }
  # if (missing(latents)){
  #   latents <- paste0("Eta_between_",seq_len(nLat))
  # }
  # if (length(latents) != nLat){
  #   stop("Length of 'latents' is not equal to number of latent variables in model.")
  # }
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLat))
  }
  if (length(latents) != nLat){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }
  
  # Generate model object:
  model <- generate_psychonetrics(model = "dlvm1", 
                                  types = list(
                                    within_latent = within_latent, between_latent = between_latent,
                                    within_residual = within_residual, between_residual = between_residual
                                  ),
                                  sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  identification=identification, verbose = verbose)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  nAllVar <- length(allVars)
  
  # Add number of observations:
  model@sample@nobs <-  
    nAllVar * (nAllVar+1) / 2 * nGroup + # Covariances per group
    nAllVar * nGroup # Means per group
  
  # Model matrices:
  modMatrices <- list()
  
  # Expected means:
  expMeans <- lapply(model@sample@means, function(x)tapply(x,timePerVar,mean,na.rm=TRUE))
  
  # Fix nu
  modMatrices$nu <- matrixsetup_mu(nu,nNode = nVar, nGroup = nGroup, labels = varnames,equal = "nu" %in% equal,
                                    expmeans = expMeans, sampletable = sampleStats, name = "nu")
  
  # Fix between-subject means:
  modMatrices$mu_eta <- matrixsetup_mu(mu_eta,nNode = nLat, nGroup = nGroup, labels = latents, equal = "mu_eta" %in% equal,
                                       expmeans = lapply(1:nGroup,function(x)rep(0, nLat)), sampletable = sampleStats, name = "mu_eta")
  
  # FIXME: simple start values for now...
  
  # T=1 cov structure:
  firstVars <- apply(vars,1,function(x)na.omit(x)[1])
  secondVars <- apply(vars,1,function(x)na.omit(x)[2])
  firstSigma0 <- lapply(sampleStats@covs,function(x)spectralshift(x[firstVars,firstVars]))
  firstSigma1 <- lapply(sampleStats@covs,function(x)spectralshift(x[secondVars,firstVars]))
  
  # Setup lambda:
  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov=firstSigma0, nGroup = nGroup, observednames = varnames, latentnames = latents,
                                           sampletable = sampleStats, name = "lambda")
  
  # If beta = 0, these sort of estimate the within and between subject covs:
  prior_bet_cov <- lapply(firstSigma1,function(x)spectralshift(0.5*(x+t(x))))
  prior_wit_cov <- lapply(seq_along(firstSigma1),function(i)spectralshift(firstSigma0[[i]] - prior_bet_cov[[i]]))
  
  
  # # Setup lambda_within:
  # modMatrices$lambda_within <- matrixsetup_lambda(lambda_within, expcov=prior_wit_cov, nGroup = nGroup, observednames = sampleStats@variables$label, latentnames = latents,
  #                                          sampletable = sampleStats, name = "lambda_within")
  
  
  # Quick and dirty sigma_zeta_within estimate:
  prior_sig_zeta_within <- lapply(seq_along(firstSigma1),function(i){
    # Let's take a pseudoinverse:
    curLam <- matrix(as.vector(modMatrices$lambda$start[,,i,drop=FALSE]),nVar,nLat)
    
    inv <- corpcor::pseudoinverse(as.matrix(kronecker(curLam,curLam)))
    
    # And obtain psi estimate:
    matrix(inv %*% as.vector(prior_wit_cov[[i]])/2,nLat,nLat)
  })
  
  
  # Setup latent varcov:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma = sigma_zeta_within,lowertri = lowertri_zeta_within,omega = omega_zeta_within,delta = delta_zeta_within,kappa = kappa_zeta_within,
                                       type = within_latent,
                                       name= "zeta_within",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nLat,
                                       expCov = prior_sig_zeta_within,
                                       nGroup = nGroup,
                                       labels = latents
                   ))
  
  # Setup Beta:
  modMatrices$beta <- matrixsetup_beta(beta, 
                                       name = "beta",
                                       nNode = nLat, 
                                       nGroup = nGroup, 
                                       labels = latents,
                                       equal = "beta" %in% equal, sampletable = sampleStats)
  
  
  # Setup residuals:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma_epsilon_within,lowertri_epsilon_within,omega_epsilon_within,delta_epsilon_within,kappa_epsilon_within,
                                       type = within_residual,
                                       name= "epsilon_within",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nVar,
                                       expCov = lapply(1:nGroup,function(x)diag(0.5,nVar)),
                                       nGroup = nGroup,
                                       labels = varnames
                   ))
  
  # Between-case effects:
  # Setup lambda_between:
  # modMatrices$lambda_between <- matrixsetup_lambda(lambda_between, expcov=prior_bet_cov, nGroup = nGroup, observednames = sampleStats@variables$label, latentnames = latents,
  #                                                 sampletable = sampleStats, name = "lambda_between")
  # 
  
  
  # Quick and dirty sigma_zeta_between estimate:
  prior_sig_zeta_between <- lapply(seq_along(firstSigma1),function(i){
    # Let's take a pseudoinverse:
    curLam <- matrix(as.vector(modMatrices$lambda$start[,,i,drop=FALSE]),nVar,nLat)
    
    inv <- corpcor::pseudoinverse(as.matrix(kronecker(curLam,curLam)))
    
    # And obtain psi estimate:
    matrix(inv %*% as.vector(prior_bet_cov[[i]])/2,nLat,nLat)
  })
  
  
  # Setup latent varcov:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between,
                                       type = between_latent,
                                       name= "zeta_between",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nLat,
                                       expCov = prior_sig_zeta_between,
                                       nGroup = nGroup,
                                       labels = latents
                   ))
  
  # Setup latent residual:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma_epsilon_between,lowertri_epsilon_between,omega_epsilon_between,delta_epsilon_between,kappa_epsilon_between,
                                       type = between_residual,
                                       name= "epsilon_between",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nVar,
                                       expCov = lapply(1:nGroup,function(x)diag(0.5,nVar)),
                                       nGroup = nGroup,
                                       labels = varnames
                   ))
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  model@extramatrices <- list(
    # Entire duplication matrix needed for likelihood:
    D = psychonetrics::duplicationMatrix(nAllVar),
    
    # Duplication matrices:
    D_y = psychonetrics::duplicationMatrix(nVar),
    D_eta = psychonetrics::duplicationMatrix(nLat),
    # D_within = psychonetrics::duplicationMatrix(nLat),
    # D_between = psychonetrics::duplicationMatrix(nLat),
    
    # Strict duplication matrices:
    Dstar_y = psychonetrics::duplicationMatrix(nVar,diag = FALSE),
    Dstar_eta = psychonetrics::duplicationMatrix(nLat,diag = FALSE),
    # Dstar_within = psychonetrics::duplicationMatrix(nLat,diag = FALSE),
    # Dstar_between = psychonetrics::duplicationMatrix(nLat,diag = FALSE),  
    
    # Identity matrices:
    I_y = as(diag(nVar),"dgCMatrix"),
    I_eta = as(diag(nLat),"dgCMatrix"),
    # I_within = Diagonal(nLat),
    # I_between = Diagonal(nLat),
    
    # Diagonalization matrices:
    A_y = psychonetrics::diagonalizationMatrix(nVar),
    A_eta = psychonetrics::diagonalizationMatrix(nLat),
    # A_within = psychonetrics::diagonalizationMatrix(nLat),
    # A_between = psychonetrics::diagonalizationMatrix(nLat),
    
    # Commutation matrices:
    C_y_y = as(lavaan::lav_matrix_commutation(nVar, nVar),"dgCMatrix"),
    C_y_eta = as(lavaan::lav_matrix_commutation(nVar, nLat),"dgCMatrix"),
    C_eta_eta = as(lavaan::lav_matrix_commutation(nLat, nLat),"dgCMatrix"),
    
    # 
    # C_y_within = as(lavaan::lav_matrix_commutation(nVar, nLat),"sparseMatrix"),
    # C_within_within = as(lavaan::lav_matrix_commutation(nLat, nLat),"sparseMatrix"),
    # C_y_between = as(lavaan::lav_matrix_commutation(nVar, nLat),"sparseMatrix"),
    # C_between_between = as(lavaan::lav_matrix_commutation(nLat, nLat),"sparseMatrix"),
    # 
    # Elimination matrices:
    L_y = psychonetrics::eliminationMatrix(nVar),
    L_eta = psychonetrics::eliminationMatrix(nLat),
    
    # L_within = psychonetrics::eliminationMatrix(nLat),
    # L_between = psychonetrics::eliminationMatrix(nLat),
    
    design = design
  )
  
  
  # Come up with P...
  # Dummy matrix to contain indices:
  # Dummy matrices with indices:
  muDummy <- matrix(rep(1:nVar,nTime))
  sigDummy <- matrix(0,nVar,nVar)
  sigDummy[lower.tri(sigDummy,diag=TRUE)] <- max(muDummy) + seq_len(nVar*(nVar+1)/2)
  sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]
  
  U <- list(sigDummy)
  # Now make all lag-k blocks...
  # Form blocks:
  for (i in 1:(nTime-1)){
    U[[length(U) + 1]] <- matrix(max(unlist(U)) + seq_len(nVar^2), nVar, nVar)
  }
  
  allSigmas <- blockToeplitz(U)
  
  # Total number:
  totElements <- max(allSigmas)
  
  # Now subset with only observed:
  subMu <- muDummy[as.vector(design==1),,drop=FALSE]
  subSigmas <- allSigmas[as.vector(design==1),as.vector(design==1)]
  
  inds <- c(as.vector(subMu),subSigmas[lower.tri(subSigmas,diag=TRUE)])
  
  # P matrix:
  # P <- bdiag(Diagonal(nVar*2),sparseMatrix(j=seq_along(inds),i=inds))
  distVec <-  c(as.vector(subMu),subSigmas[lower.tri(subSigmas,diag=TRUE)])
  nTotal <- length(distVec)
  distVecrawts <- seq_along(distVec)[distVec!=0]
  distVec <- distVec[distVec!=0]
  # Now I can make the matrix:
  
  model@extramatrices$P <- sparseMatrix(
    i = distVecrawts, j = distVec, dims = c(nTotal, totElements)
  )
  
  model@extramatrices$P <- as( model@extramatrices$P, "dgCMatrix")
  # model@extramatrices$P <- sparseMatrix(j=seq_along(inds),i=order(inds))
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){
    
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
                                                baseline_saturated = FALSE,
                                                sampleStats = sampleStats)
    
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
                                                 baseline_saturated = FALSE,
                                                 sampleStats = sampleStats)
    
    # if not FIML, Treat as computed:
    if (estimator != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)      
    }
  }
  
  # Identify model:
  if (identify){
    model <- identify(model)
  }
  
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }
  
  # Return model:
  return(model)
}
