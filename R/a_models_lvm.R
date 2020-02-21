# Latent variable model (lvm)
lvm <- function(
  data, # Dataset
  lambda, # only required non-missing matrix
  latent = c("cov","chol","prec","ggm"), # Maybe add cor at some point, but not now
  residual = c("cov","chol","prec","ggm"),
  
  # Latent matrices:
  sigma_zeta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa_zeta = "full", # Precision
  omega_zeta = "full", # Partial correlations
  lowertri_zeta = "full", # Cholesky
  delta_zeta = "full", # Used for both ggm and pcor
  
  # Residual matrices:
  sigma_epsilon = "empty", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa_epsilon = "empty", # Precision
  omega_epsilon = "empty", # Partial correlations
  lowertri_epsilon = "empty", # Cholesky
  delta_epsilon = "empty", # Used for both ggm and pcor
  
  # Beta:
  beta = "empty",
  
  # Mean structure:
  nu,
  nu_eta,
  
  # Identification:
  identify = TRUE,
  identification = c("loadings","variance"),
  
  # Rest:
  vars, # character indicating the variables Extracted if missing from data - group variable
  latents, # Name of latent varianles
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
  covtype = c("choose","ML","UB"),
  standardize = c("none","z","quantile"),
  sampleStats,
  verbose=TRUE
){
  rawts = FALSE
  if (rawts){
    warning("'rawts' is only included for testing purposes. Please do not use!")
  }
  
  # Type:
  latent <- match.arg(latent)
  residual <- match.arg(residual)
  
  # Identification:
  identification <- match.arg(identification)

  # WLS weights:
  if (missing(WLS.V)){
    WLS.V <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                    switch(estimator,
                           "WLS" = "full",
                           "ULS" = "identity",
                           "DWLS" = "diag"
                    ))
  }
  
  # Obtain sample stats:
  if (missing(sampleStats)){
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
                               covtype=covtype,
                               weightsmatrix = WLS.V,
                               verbose=verbose,
                               standardize=standardize)
  }


  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "lvm",sample = sampleStats,computed = FALSE, 
                                  equal = equal,identification=identification,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(latent = latent, residual = residual))
  
  # Submodel:
  latentCov <- latent %in% c("cov","chol")
  residCov <- residual %in% c("cov", "chol")
  
  if (latentCov & residCov){
    model@submodel <- "sem"
  } else if (latentCov & !residCov){
    model@submodel <- "rnm"
  } else if (!latentCov & residCov){
    model@submodel <- "lnm"
  } else if (!latentCov & !residCov){
    model@submodel <- "lrnm"
  }
  
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
  
  # Number of latents:
  nLatent <- ncol(lambda)
  
  
  # If latents is not provided, make it:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }
  
  # Model matrices:
  modMatrices <- list()
  
  # Fix nu
  modMatrices$nu <- matrixsetup_mu(nu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "nu" %in% equal,
                       expmeans = model@sample@means, sampletable = sampleStats, name = "nu")
  
  # Fix nu_eta
  modMatrices$nu_eta <- matrixsetup_mu(nu_eta,nNode = nLatent,nGroup = nGroup,labels = latents,equal = "nu_eta" %in% equal,
                                    expmeans = lapply(seq_len(nGroup),function(x)rep(0,nLatent)), sampletable = sampleStats, name = "nu_eta")
  
   # Setup lambda:
  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov=model@sample@covs, nGroup = nGroup, observednames = sampleStats@variables$label, latentnames = latents, sampletable = sampleStats)
  
  # Setup beta:
  modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup, labels = latents, sampletable = sampleStats)
  
  # Compute the expected latent and residual cov matrices:
  expLatSigma <- lapply(1:nGroup,function(x)matrix(0,nLatent,nLatent))
  expResidSigma <- lapply(1:nGroup,function(x)matrix(0,nNode,nNode))
  
  # For each group:
  for (g in 1:nGroup){
    # Current cov estimate:
    curcov <- as.matrix(sampleStats@covs[[g]])
 
    # Cur loadings:
    curLambda <- modMatrices$lambda$start[,,g]
    if (!is.matrix(curLambda)){
      curLambda <- as.matrix(curLambda)
    }
    
    # Residual variances, let's start by putting the vars on 1/4 times the observed variances:
    Theta <- diag(diag(curcov)/4)
    
    # Check if this is positive definite:
    ev <- eigen(curcov - Theta)$values
    
    # Shrink until it is positive definite:
    loop <- 0
    repeat{
      ev <- eigen(curcov - Theta)$values
      if (loop == 100){
        # give up...
        
        Theta <- diag(nrow(Theta))
        break
      }
      if (all(ev>0)){
        break
      }
      Theta <- Theta * 0.9
      loop <- loop + 1
    }
    
    # Expected residual sigma:
    expResidSigma[[g]] <- Theta
    
    # This means that the factor-part is expected to be:
    factorPart <- curcov - Theta
    
    # Let's take a pseudoinverse:
    inv <- corpcor::pseudoinverse(kronecker(curLambda,curLambda))
    
    # And obtain psi estimate:
    expLatSigma[[g]] <- matrix(inv %*% as.vector(factorPart),nLatent,nLatent)
  }
  
  # Latent varcov:
  if (latent == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta, 
                                                name = "sigma_zeta",
                                           expcov=expLatSigma,
                                           nNode = nLatent, 
                                           nGroup = nGroup, 
                                           labels = latents,
                                           equal = "sigma_zeta" %in% equal, sampletable = sampleStats,
                                           beta = modMatrices$beta[[1]])    
  } else if (latent == "chol"){
    modMatrices$lowertri_zeta <- matrixsetup_lowertri(lowertri_zeta, 
                                                      name = "lowertri_zeta",
                                                      expcov=expLatSigma,
                                                      nNode = nLatent, 
                                                      nGroup = nGroup, 
                                                      labels = latents,
                                                 equal = "lowertri_zeta" %in% equal, sampletable = sampleStats,
                                                 beta = modMatrices$beta[[1]])
  } else if (latent == "ggm"){
    # Add omega matrix:
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta, 
                                                name = "omega_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent, 
                                                nGroup = nGroup, 
                                                labels = latents,
                                           equal = "lowertri_zeta" %in% equal, sampletable = sampleStats,
                                           beta = modMatrices$beta[[1]])
    
    # Add delta matrix:
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta, 
                                                name = "delta_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent, 
                                                nGroup = nGroup, 
                                                labels = latents,
                                           equal = "delta_zeta" %in% equal, sampletable = sampleStats) 
  } else if (latent == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_zeta <- matrixsetup_kappa(kappa_zeta, 
                                                   name = "kappa_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent, 
                                                nGroup = nGroup, 
                                                labels = latents,
                                           equal = "kappa_zeta" %in% equal, sampletable = sampleStats,
                                           beta = modMatrices$beta[[1]])
  }
  
  ### Residual varcov ###
  if (residual == "cov"){
    modMatrices$sigma_epsilon <- matrixsetup_sigma(sigma_epsilon, 
                                                   name = "sigma_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label,
                                                equal = "sigma_epsilon" %in% equal, sampletable = sampleStats)    
  } else if (residual == "chol"){
    modMatrices$lowertri_epsilon <- matrixsetup_lowertri(lowertri_epsilon, 
                                                         name = "lowertri_epsilon",
                                                      expcov=expResidSigma,
                                                      nNode = nNode, 
                                                      nGroup = nGroup, 
                                                      labels = sampleStats@variables$label,
                                                      equal = "lowertri_epsilon" %in% equal, sampletable = sampleStats)
  } else if (residual == "ggm"){
    # Add omega matrix:
    modMatrices$omega_epsilon <- matrixsetup_omega(omega_epsilon, 
                                                   name = "omega_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label,
                                                equal = "omega_epsilon" %in% equal, sampletable = sampleStats)
    
    # Add delta matrix:
    modMatrices$delta_epsilon <- matrixsetup_delta(delta_epsilon, 
                                                   name = "delta_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label,
                                                equal = "delta_epsilon" %in% equal, sampletable = sampleStats) 
  } else if (residual == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_epsilon <- matrixsetup_kappa(kappa_epsilon, 
                                                   name = "kappa_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label,
                                                equal = "kappa_epsilon" %in% equal, sampletable = sampleStats)
  }
  
  
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  # Extra matrices:
  model@extramatrices <- list(
    D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    Deta = psychonetrics::duplicationMatrix(nLatent), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    L_eta = psychonetrics::eliminationMatrix(nLatent), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    Dstar_eta = psychonetrics::duplicationMatrix(nLatent,diag = FALSE), # Strict duplicaton matrix
    In = Diagonal(nNode), # Identity of dim n
    Inlatent = Diagonal(nLatent),
    C = as(lavaan::lav_matrix_commutation(nNode, nLatent),"sparseMatrix"),
    Cbeta = as(lavaan::lav_matrix_commutation(nLatent, nLatent),"sparseMatrix"),
    C_chol = as(lavaan::lav_matrix_commutation(nNode, nNode),"sparseMatrix"),
    A = psychonetrics::diagonalizationMatrix(nNode),
    Aeta = psychonetrics::diagonalizationMatrix(nLatent)
  )
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
   
    # Form baseline model:
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
                                             baseline_saturated = FALSE, sampleStats = sampleStats)
    
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
           baseline_saturated = FALSE, sampleStats = sampleStats)
    
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
  
  # Return model:
  return(model)
}