# Latent network model creator
ml_lvm <- function(
  data, # Dataset
  
  # Factor loadings:
  lambda, # May not be missing
  clusters, # Cluster variable, may not be missing
  # lambda_within, # May not be missing
  # lambda_between, # May not be missing
  
  # Type:
  within_latent = c("cov","chol","prec","ggm"), 
  within_residual = c("cov","chol","prec","ggm"), 
  between_latent = c("cov","chol","prec","ggm"), 
  between_residual = c("cov","chol","prec","ggm"), 
  
  # Beta:
  beta_within = "empty",
  beta_between = "empty",
  
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
  nu_eta,
  
  # Identification:
  identify = TRUE,
  identification = c("loadings","variance"),
  
  vars, # character indicating the variables Extracted if missing from data - group variable
  
  latents,
  
  # The rest:
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  # fitfunctions, # Leave empty
  estimator = c("FIML","MUML"),
  optimizer,
  storedata = FALSE,
  verbose = FALSE,
  standardize = c("none","z","quantile"),
  sampleStats
){

  # CRAN Check workarounds (sorry):
  . <- NULL
  variable <- NULL
  value <- NULL
  
  # Check input:
  # Match args:
  standardize <- match.arg(standardize)
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  within_residual <- match.arg(within_residual)
  between_residual <- match.arg(between_residual)
  identification <- match.arg(identification)
  estimator <- match.arg(estimator)
  if (estimator != "FIML") stop("Only FIML supported.")
  
  # Check clusters:
  if (missing(clusters)){
    stop("'clusters' may not be missing (use lvm family instead).")
  }
  
  # Check if data is data frame:
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  
  # Check vars:
  if (missing(vars)){
    if (is.null(names(data))){
      stop("Dataset contains no column names.")
    }
    vars <- names(data)
    vars <- vars[vars!=clusters]
  } else {
    if (is.null(names(data))){
      stop("Dataset contains no column names.")
    }
  }
  
  # Number of variables:
  nVar <- length(vars)
  
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
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  # Number of clusters:
  if (!clusters %in% names(data)){
    stop("'clusters' argument does not correspond to column name of 'data'")
  }
  
  # Remove data with NA cluster:
  if (any(is.na(data[[clusters]]))){
    warning("Rows with NA cluster removed.")
    data <- data[!is.na(data[[clusters]]),]
  }
  
  # Number of clusters:
  nCluster <- length(unique(data[[clusters]]))
  
  # Add a column with ID in cluster:
  data[['CLUSTERID']] <- unlist(tapply(data[[clusters]],data[[clusters]],seq_along))
  
  # Add group:
  if (missing(groups)){
    groups <- "GROUPID"
    data[[groups]] <- "fullsample"
  }
  
  # Max in cluster:
  maxInCluster <- max(data[['CLUSTERID']])
  
  
  # Standardize the data:
  if (standardize == "z"){
    for (v in seq_along(vars)){
      data[,vars[v]] <- (data[,vars[v]] - mean(data[,vars[v]],na.rm=TRUE)) / sd(data[,vars[v]],na.rm=TRUE)
    }
  } else if (standardize == "quantile"){
    for (v in seq_along(vars)){
      data[,vars[v]] <- quantiletransform(data[,vars[v]])
    }
  }
  
  ### For start values, compute pairwise covs of sample means and deviations from means ###
  
  # Mean dataset:
  data_means <- data %>% group_by(!!as.name(clusters),!!as.name(groups)) %>% summarize_at(.vars=vars,funs(mean(.,na.rm=TRUE)))
  
  # Within person dataset:
  data_centered <- data %>% group_by(!!as.name(clusters),!!as.name(groups)) %>% 
    mutate_at(.vars=vars,funs(scale(.,scale = FALSE)))
  
  
  start_between_covs <- list()
  start_within_covs <- list()
  start_covs <- list()
  
  allGroups <- unique(data[[groups]])
  for (i in seq_along(allGroups)){
    start_between_covs[[i]] <- cov(data_means[data_means[[groups]] == allGroups[i],vars],use="pairwise.complete.obs")
    start_within_covs[[i]] <- cov(data_centered[data_centered[[groups]] == allGroups[i],vars],use="pairwise.complete.obs")
    start_covs[[i]] <- cov(data[,vars],use="pairwise.complete.obs")
  }
  
  ###
  
  # To long format:
  datalong <- tidyr::gather(data,variable,value,vars)
  
  # Transform data to wide format:
  datawide <- tidyr::pivot_wider(datalong, id_cols = c(clusters,groups),  values_from = "value", names_from = c("variable","CLUSTERID"))
  
  # Now make a design matrix:
  rowVars <- vars
  colVars <- as.character(seq(maxInCluster))
  design <- matrix(NA, length(rowVars), length(colVars))
  for (i in seq_along(rowVars)){
    for (j in seq_along(colVars)){
      varName <- paste0("^",rowVars[i],"_",colVars[j],"$")
      whichVar <- which(grepl(varName, names(datawide)))
      if (length(whichVar) == 1){
        design[i,j] <- paste0(rowVars[i],"_",colVars[j])
      }
    }
  }
  
  
  # List all variables to use, in order:
  allVars <- na.omit(as.vector(design))
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    sampleStats <- samplestats(data = datawide, 
                               vars = allVars, 
                               groups = groups,
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               verbose=verbose)
  }
  
  
  
  # designPattern matrix:
  designPattern <- as(1*(!is.na(design)),"matrix")
  
  # cases per var:
  casesPerVar <- as.vector(designPattern * row(designPattern))
  casesPerVar <- casesPerVar[casesPerVar!=0]
  
  
  # row names:
  if (is.null(rownames(design))){
    rownames(design) <- paste0("V",seq_len(nrow(design)))
  }
  varnames <- rownames(design)
  
  # col names:
  if (is.null(colnames(design))){
    colnames(design) <- paste0("C",seq_len(ncol(design)))
  }
  casenames <- colnames(design)
  
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
  model <- generate_psychonetrics(model = "ml_lvm", 
                                  types = list(
                                    within_latent = within_latent, between_latent = between_latent,
                                    within_residual = within_residual, between_residual = between_residual
                                  ),
                                  sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  identification=identification, verbose=verbose)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Number of total vars:
  nAllVar <- length(allVars)
  
  # Add number of observations:
  model@sample@nobs <-  
    nAllVar * (nAllVar+1) / 2 * nGroup + # Covariances per group
    nAllVar * nGroup # Means per group
  
  # Model matrices:
  modMatrices <- list()
  
  # Expected means:
  expMeans <- lapply(model@sample@means, function(x)tapply(x,casesPerVar,mean,na.rm=TRUE))
  
  # Fix nu
  modMatrices$nu <- matrixsetup_mu(nu,nNode = nVar, nGroup = nGroup, labels = varnames,equal = "nu" %in% equal,
                                   expmeans = expMeans, sampletable = sampleStats, name = "nu")
  
  # Fix between-subject means:
  modMatrices$nu_eta <- matrixsetup_mu(nu_eta,nNode = nLat, nGroup = nGroup, labels = latents, equal = "nu_eta" %in% equal,
                                       expmeans = lapply(1:nGroup,function(x)rep(0, nLat)), sampletable = sampleStats, name = "nu_eta")
  
  
  # Setup lambda:
  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov=start_covs, nGroup = nGroup, 
                                           observednames = varnames, latentnames = latents,
                                           sampletable = sampleStats, name = "lambda")
  
  #### Within-person model ####
  # Setup beta:
  modMatrices$beta_within <- matrixsetup_beta(beta_within,name = "beta_within", nNode = nLat, nGroup = nGroup, labels = latents, sampletable = sampleStats)
  
  
  ## Prior guesses for latent and residual variances:
  priorsWithin <- expected_latent_residual_covs(start_within_covs, modMatrices$lambda$start)
  
  # Setup latent varcov:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma = sigma_zeta_within,lowertri = lowertri_zeta_within,omega = omega_zeta_within,delta = delta_zeta_within,kappa = kappa_zeta_within,
                                       type = within_latent,
                                       name= "zeta_within",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nLat,
                                       expCov = priorsWithin$latent,
                                       nGroup = nGroup,
                                       labels = latents
                   ))
  
  
  
  # Setup residuals:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma_epsilon_within,lowertri_epsilon_within,omega_epsilon_within,delta_epsilon_within,kappa_epsilon_within,
                                       type = within_residual,
                                       name= "epsilon_within",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nVar,
                                       expCov = priorsWithin$residual,
                                       nGroup = nGroup,
                                       labels = varnames
                   ))
  
  
  #### Between-person model ####
  # Setup beta:
  modMatrices$beta_between <- matrixsetup_beta(beta_between,name = "beta_between", nNode = nLat, nGroup = nGroup, labels = latents, sampletable = sampleStats)
  
  ## Prior guesses for latent and residual variances:
  priorsBetween <- expected_latent_residual_covs(start_between_covs, modMatrices$lambda$start)
  
  # Setup latent varcov:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between,
                                       type = between_latent,
                                       name= "zeta_between",
                                       sampleStats= sampleStats,
                                       equal = equal,
                                       nNode = nLat,
                                       expCov = priorsBetween$latent,
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
                                       expCov = priorsBetween$residual,
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
    
    designPattern = designPattern
  )

  # Come up with P...
  # Dummy matrix to contain indices:
  # Dummy matrices with indices:
  muDummy <- matrix(rep(1:nVar,maxInCluster))
  sigDummy <- matrix(0,nVar,nVar)
  sigDummy[lower.tri(sigDummy,diag=TRUE)] <- max(muDummy) + seq_len(nVar*(nVar+1)/2)
  sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]
  
  U <- list(sigDummy)
  # Now make all lag-k blocks...
  # Form blocks:
  if (maxInCluster > 1){
    U <- c(U,lapply(seq_len(maxInCluster-1),function(x)matrix(max(sigDummy) + seq_len(nVar^2), nVar, nVar)))
  }
  
  
  allSigmas <- blockToeplitz(U)
  
  # Total number:
  totElements <- max(allSigmas)
  
  # Now subset with only observed:
  subMu <- muDummy[as.vector(designPattern==1),,drop=FALSE]
  subSigmas <- allSigmas[as.vector(designPattern==1),as.vector(designPattern==1)]
  
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
  
  model@extramatrices$P <- as(model@extramatrices$P, "dgCMatrix")
  # model@extramatrices$P <- sparseMatrix(j=seq_along(inds),i=order(inds))
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){
    
    # Form baseline model:
    # model@baseline_saturated$baseline <- varcov(data,
    #                                             type = "chol",
    #                                             lowertri = "empty",
    #                                             vars = allVars,
    #                                             groups = groups,
    #                                             covs = covs,
    #                                             means = means,
    #                                             nobs = nobs,
    #                                             missing = missing,
    #                                             equal = equal,
    #                                             estimator = estimator,
    #                                             baseline_saturated = FALSE,
    #                                             sampleStats = sampleStats)
    # within_latent = c("cov","chol","prec","ggm"), 
    # within_residual = c("cov","chol","prec","ggm"), 
    # between_latent = c("cov","chol","prec","ggm"), 
    # between_residual = c("cov","chol","prec","ggm"), 
    # 
    I <- diag(nVar)
    O <- matrix(0, nVar, nVar)
    model@baseline_saturated$baseline <- ml_lvm(data,
                                                within_latent = "chol",
                                                lowertri_zeta_within = "empty",
                                                between_latent = "chol",
                                                lowertri_zeta_between = "empty",
                                                lowertri_epsilon_within = O,
                                                lowertri_epsilon_between = O,
                                                lambda = I,
                                                vars = vars,
                                                clusters = clusters,
                                                groups = groups,
                                                equal = equal,
                                                estimator = estimator,
                                                baseline_saturated = FALSE,
                                                sampleStats = sampleStats)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- ml_lvm(data,
                                                 within_latent = "chol",
                                                 sigma_zeta_within = "full",
                                                 between_latent = "chol",
                                                 sigma_zeta_between = "full",
                                                 lowertri_epsilon_within = O,
                                                 lowertri_epsilon_between = O,
                                                 lambda = I,
                                                 vars = vars,
                                                 clusters = clusters,
                                                 groups = groups,
                                                 equal = equal,
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats = sampleStats)
    
    
    # model@baseline_saturated$saturated <- varcov(data,
    #                                              type = "chol", 
    #                                              lowertri = "full", 
    #                                              vars = allVars,
    #                                              groups = groups,
    #                                              covs = covs,
    #                                              means = means,
    #                                              nobs = nobs,
    #                                              missing = missing,
    #                                              equal = equal,
    #                                              estimator = estimator,
    #                                              baseline_saturated = FALSE,
    #                                              sampleStats = sampleStats)
    
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
