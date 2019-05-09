# Latent network model creator
var1 <- function(
  data, # Dataset
  
  # Type:
  contemporaneous = c("cov","chol","prec","ggm"), 
  
  # Temporal effects:
  beta = "full",
  
  # Contemporaneous effects:
  omega_zeta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_zeta = "full", # If missing, just full for both groups or equal
  kappa_zeta = "full",
  sigma_zeta = "full",
  lowertri_zeta = "full",
  
  # The rest:
  mu,
  beepvar,
  dayvar,
  idvar,
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  # fitfunctions, # Leave empty
  estimator = "ML",
  optimizer = "default", # ucminf
  storedata = FALSE,
  sampleStats
){
  contemporaneous <- match.arg(contemporaneous)
  
  # FIXME: Not sure why needed...
  if (missing(vars)) vars2 <- NULL else vars2 <- vars
  if (missing(idvar)) idvar <- NULL
  if (missing(dayvar)) dayvar <- NULL
  if (missing(beepvar)) beepvar <- NULL
  if (missing(groups)) groups <- NULL
  
  # If data is missing with rawts, stop:
  if (!missing(data)){
    data <- as.data.frame(data)
    if (is.null(names(data))){
      names(data) <- paste0("V",seq_len(ncol(data)))
    }
  }
  
  # If data is not missing, make augmented data:
  data <- tsData(data, vars = vars2, beepvar = beepvar, dayvar = dayvar, idvar = idvar, groupvar = groups)
  
  # Extract var names:
  if (is.null(groups)){
    vars <- colnames(data)  
  } else {
    vars <- colnames(data)[colnames(data)!=groups]
  }
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    sampleStats <- samplestats(data = data, 
                               vars = vars, 
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing  = ifelse(estimator == "FIML","pairwise",missing),
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               weightsmatrix = ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                                                      switch(estimator,
                                                        "WLS" = "full",
                                                        "ULS" = "identity",
                                                        "DWLS" = "diag"
                                                      )))    
  }

  
  # Check if number of variables is an even number:
  if ( nrow(sampleStats@variables) %% 2 != 0){
    stop("Number of variables is not an even number: variance-covariance matrix cannot be a Toeplitz matrix. ")
  }
  
  
  # Check some things:
  nNode <- nrow(sampleStats@variables) / 2
  
  
  # Generate model object:
  model <- generate_psychonetrics(model = "var1",submodel = 
                                    switch(contemporaneous,
                                           "prec" = "gvar",
                                           "ggm" = "gvar",
                                           "chol" = "var",
                                           "cov" = "var"
                                    ),types = list(zeta = contemporaneous),
                                  sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # FIXME: Keep this the same for now for rawts = TRUE
  nVar <- nNode * 2
  # Add number of observations:
  model@sample@nobs <-  
    nVar * (nVar+1) / 2 * nGroup + # Covariances per group
    nVar * nGroup # Means per group
  
  # Model matrices:
  modMatrices <- list()
  
  # Fix mu
  modMatrices$mu <- matrixsetup_mu(mu,nNode = nVar, nGroup = nGroup, labels = sampleStats@variables$label,equal = "mu" %in% equal,
                                   expmeans = model@sample@means, sampletable = sampleStats, name = "mu")
  
  shiftCovs <- lapply(sampleStats@covs,spectralshift)
  
  # Exogenous block covariances:
  exoCovs <- lapply(shiftCovs,function(x)spectralshift(x[1:nNode,1:nNode]))
  
  # Fix exo cholesky:
  modMatrices$exo_cholesky <- matrixsetup_lowertri("full", 
                                                   name = "exo_cholesky",
                                                   expcov=exoCovs,
                                                   nNode = nNode, 
                                                   nGroup = nGroup, 
                                                   labels = sampleStats@variables$label[1:nNode],
                                                   equal = "exo_cholesky" %in% equal, sampletable = sampleStats)
  
  
  # S1 and S0 estimates:
  S0est <- lapply(shiftCovs,function(x)spectralshift(x[nNode + (1:nNode),nNode + (1:nNode)]))
  S1est <- lapply(shiftCovs,function(x)x[nNode + (1:nNode),1:nNode])
  S0inv <- lapply(S0est,solve_symmetric)
  
  # Prior estimate for beta:
  betaEst <- lapply(1:nGroup, function(g) as.matrix(S1est[[g]] %*% S0inv[[g]]))
  
  modMatrices$beta <- matrixsetup_beta(beta, 
                                       name = "beta",
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = sampleStats@variables$label[nNode + (1:nNode)],
                                       equal = "beta" %in% equal, sampletable = sampleStats, start = betaEst,
                                       onlyStartSign = FALSE)
  
  
  # A prior guess for the contemporaneous covariances is (Schur complement):
  # contCovEst <- lapply(1:nGroup, function(g) spectralshift(S0est[[g]] - t(S1est[[g]]) %*% S0inv[[g]] %*% S1est[[g]]))
  contCovEst <- lapply(1:nGroup, function(g) spectralshift(S0est[[g]] - S1est[[g]] %*% S0inv[[g]] %*% t(S1est[[g]])))
  # contCovEst <- lapply(1:nGroup, function(g) spectralshift(exoCovs[[g]] - t(S1est[[g]]) %*% S0inv[[g]] %*% S1est[[g]]))
  
  # Fill in:
  if (contemporaneous == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta, name = "sigma_zeta",
                                                expcov=contCovEst,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label[-(1:nNode)],
                                                equal = "sigma_zeta" %in% equal, sampletable = sampleStats)    
  } else if (contemporaneous == "chol"){
    modMatrices$lowertri_zeta <- matrixsetup_lowertri(lowertri_zeta,  name = "lowertri_zeta",
                                                      expcov=contCovEst,
                                                      nNode = nNode, 
                                                      nGroup = nGroup, 
                                                      labels = sampleStats@variables$label[-(1:nNode)],
                                                      equal = "lowertri_zeta" %in% equal, sampletable = sampleStats)
  } else if (contemporaneous == "ggm"){
    # Add omega matrix:
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta, name = "omega_zeta",
                                                expcov=contCovEst,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label[-(1:nNode)],
                                                equal = "omega_zeta" %in% equal, sampletable = sampleStats,
                                                onlyStartSign = FALSE)
    
    # Add delta matrix:
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta, name = "delta_zeta",
                                                expcov=contCovEst,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label[-(1:nNode)],
                                                equal = "delta_zeta" %in% equal, sampletable = sampleStats,
                                                onlyStartSign = FALSE) 
  } else if (contemporaneous == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_zeta <- matrixsetup_kappa(kappa_zeta, name = "kappa_zeta",
                                                expcov=contCovEst,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = sampleStats@variables$label[-(1:nNode)],
                                                equal = "kappa_zeta" %in% equal, sampletable = sampleStats)
  }
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  model@extramatrices <- list(
    D =  psychonetrics::duplicationMatrix(nNode*2), # Toeplitz matrix D 
    D2 = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    In = Diagonal(nNode), # Identity of dim n
    In2 = Diagonal(nNode), # Identity of dim n^2
    A = psychonetrics::diagonalizationMatrix(nNode),
    C = as(lavaan::lav_matrix_commutation(nNode,nNode),"sparseMatrix")
    # P=P # Permutation matrix
  )
  
  
    # Come up with P...
    # Dummy matrix to contain indices:
    dummySigma <- matrix(0,nNode*2,nNode*2)
    smallMat <- matrix(0,nNode,nNode)
    dummySigma[1:nNode,1:nNode][lower.tri(smallMat,diag=TRUE)] <- seq_len(nNode*(nNode+1)/2)
    dummySigma[nNode + (1:nNode),nNode + (1:nNode)][lower.tri(smallMat,diag=TRUE)] <- max(dummySigma) + seq_len(nNode*(nNode+1)/2)
    dummySigma[nNode + (1:nNode),1:nNode] <- max(dummySigma) + seq_len(nNode^2)
    inds <- dummySigma[lower.tri(dummySigma,diag=TRUE)]
    
    # P matrix:
    # P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=inds))
    model@extramatrices$P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=order(inds)))

  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Baseline GGM should be block matrix:
    basGGM <- diag(nNode*2)
    basGGM[1:nNode,1:nNode] <- 1
    
    model@baseline_saturated$baseline <- cholesky(data = data,
                                                  lowertri = basGGM,
                                                  vars = vars,
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
    model@baseline_saturated$saturated <- cholesky(data = data, 
                                                   lowertri = "full", 
                                                   vars = vars,
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
  
  
  # Return model:
  return(model)
}
