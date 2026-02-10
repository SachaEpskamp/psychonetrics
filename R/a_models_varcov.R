# Latent network model creator
varcov <- function(
  data, # Dataset
  type = c("cov","chol","prec","ggm","cor"),
  sigma = "full", 
  kappa = "full", # Precision
  # rho = "full", # Correlations
  omega = "full", # Partial correlations
  lowertri = "full", # Cholesky
  delta = "diag", # Used for both ggm and pcor
  rho = "full", # Used for cor
  SD = "full", # Used for cor
  mu,
  tau,
  vars, # character indicating the variables Extracted if missing from data - group variable
  ordered = character(0), # character indicating the variables that are ordinal
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "auto",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer,
  storedata = FALSE,
  WLS.W,
  sampleStats, # Leave to missing
  meanstructure, # Defaults to TRUE if data is used or means is used, FALSE otherwie
  corinput, # Defaults to TRUE if the input is detected to consist of correlation matrix/matrices, FALSE otherwise
  verbose = FALSE,
  covtype = c("choose","ML","UB"),
  standardize = c("none","z","quantile"),
  fullFIML=FALSE,
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = 0,  # Penalty strength (used when estimator = "PML")
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
){
  rawts = FALSE
  if (rawts){
    warning("'rawts' is only included for testing purposes. Please do not use!")
  }
  
  # Reset ordered if needed:
  if (identical(ordered, FALSE)){
    ordered <- character(0)
  }
  
  # WLSMV is a synonym for DWLS (DWLS estimation + scaled test statistic):
  if (estimator == "WLSMV"){
    estimator <- "DWLS"
  }

  if (estimator == "default"){
    if (length(ordered) > 0){
      estimator <- "DWLS"
    } else {
      estimator <- "ML"
    }
  }

  # Experimental warnings:
  if (length(ordered) > 0) {
    experimentalWarning("ordinal data in varcov()")
  }

  # Check WLS for ordinal:
  if (length(ordered) > 0 & !estimator %in% c("WLS","DWLS","ULS")){
    stop("Ordinal data is only supported for WLS, DWLS and ULS estimators.")
  }
  
  # Check corinput FIXME: may be mixture
  if (length(ordered) > 0 & missing(corinput)){
    corinput <- TRUE
  }
  if (length(ordered) > 0 && (!missing(corinput) && !corinput)){
    stop("corinput must be TRUE for ordinal data.")
  }

  
  # Disable meanstructure (FIXME: Allow for both ordinal and continuous):
  if (length(ordered) > 0){
    meanstructure <- FALSE
  }
  
  # Type:
  type <- match.arg(type)
  
  # Set meanstructure:
  if (missing(meanstructure)){
    meanstructure <- (!missing(data) || !missing(means))
  }
  
  # Check FIML:
  if (!missing(data) && !meanstructure && estimator %in% c("FIML", "FIPML")){
    stop("meanstructure = FALSE is not yet supported for 'FIML'/'FIPML' estimator")
  }

  # Auto-detect missing data handling:
  if (missing == "auto") {
    if (!missing(data)) {
      if (missing(vars)) {
        check_vars <- if (!is.null(colnames(data))) colnames(data) else seq_len(ncol(data))
      } else {
        check_vars <- vars
      }
      has_missing <- any(is.na(data[, check_vars, drop = FALSE]))
      if (has_missing) {
        if (estimator == "ML") {
          estimator <- "FIML"
        } else if (estimator == "PML") {
          estimator <- "FIPML"
        } else {
          # LS variants: default to listwise (WLS weights don't support missing data for continuous)
          missing <- "listwise"
        }
      } else {
        missing <- "listwise"
      }
    } else {
      # No raw data (covs/means provided): fall back to listwise
      missing <- "listwise"
    }
  }

  # Obtain sample stats:
  if (missing(sampleStats)){
    # WLS weights:
    if (missing(WLS.W)){
      WLS.W <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                      switch(estimator,
                             "WLS" = "full",
                             "ULS" = "identity",
                             "DWLS" = "diag"
                      ))
    }
    
    sampleStats <- samplestats(data = data, 
                               vars = vars, 
                               ordered=ordered,
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing = ifelse(estimator %in% c("FIML", "FIPML"),"pairwise",missing),
                               rawts = rawts,
                               fimldata = estimator %in% c("FIML", "FIPML"),
                               storedata = storedata,
                               weightsmatrix = WLS.W,
                               meanstructure = meanstructure,
                               corinput = corinput,
                               covtype=covtype,
                               verbose=verbose,
                               standardize=standardize,
                               fullFIML=fullFIML,
                               bootstrap=bootstrap,
                               boot_sub = boot_sub,
                               boot_resample = boot_resample)
  }
 
  # Overwrite corinput:
  corinput <- sampleStats@corinput
  
  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "varcov",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(y = type),
                                  submodel = type, meanstructure = meanstructure, verbose=verbose)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  

  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  if (length(ordered) > 0){
    nThresh <- sum(sapply(model@sample@thresholds,function(x)sum(sapply(x,length))))    
  } else {
    nThresh <- 0
  }

  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode-1) / 2 * nGroup + # Covariances per group
    (!corinput) * nNode * nGroup + # Variances (ignored if correlation matrix is input)
    meanstructure * nMeans + # Means per group
    nThresh
  
  # Model matrices:
  modMatrices <- list()
  
  # # Fix mu
  # modMatrices$mu <- matrixsetup_mu(mu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "mu" %in% equal,
  #                                  expmeans = model@sample@means, sampletable = sampleStats, meanstructure = meanstructure)
  
  # Alternative, without mean structure just ignore mu:
  if (meanstructure){
    # Fix mu
    modMatrices$mu <- matrixsetup_mu(mu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "mu" %in% equal,
                       expmeans = model@sample@means, sampletable = sampleStats, meanstructure = meanstructure)
  }
  
  # Thresholds:
  if (length(ordered) > 0){
    # Ideal setup for tau is a matrix, with NA for the missing elements.
    modMatrices$tau <- matrixsetup_tau(tau, nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,
                                       equal = "tau" %in% equal, sampleThresholds = model@sample@thresholds, sampletable = sampleStats)
    
  }
  
  # fixMu(mu,nGroup,nNode,"mu" %in% equal)

  # Fix sigma
  if (type == "cov"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'cov'. Use type = 'cor' or set corinput = FALSE")
    }
    
    modMatrices$sigma <- matrixsetup_sigma(sigma, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "sigma" %in% equal, sampletable = sampleStats)    
  } else if (type == "chol"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'chol'.")
    }
    modMatrices$lowertri <- matrixsetup_lowertri(lowertri, 
                                                 expcov=model@sample@covs,
                                                 nNode = nNode, 
                                                 nGroup = nGroup, 
                                                 labels = sampleStats@variables$label,
                                                 equal = "lowertri" %in% equal, sampletable = sampleStats)
  } else if (type == "ggm"){
    # Add omega matrix:
    modMatrices$omega <- matrixsetup_omega(omega, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "omega" %in% equal, sampletable = sampleStats)
    
    if (!corinput){
      # Add delta matrix (ingored if corinput == TRUE):
      modMatrices$delta <- matrixsetup_delta(delta, 
                                             expcov=model@sample@covs,
                                             nNode = nNode, 
                                             nGroup = nGroup, 
                                             labels = sampleStats@variables$label,
                                             equal = "delta" %in% equal, sampletable = sampleStats,
                                             omegaStart =  modMatrices$omega$start)       
    }
    
  } else if (type == "prec"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'prec'. Use type = 'ggm' or set corinput = FALSE")
    }
    
    # Add omega matrix:
    modMatrices$kappa <- matrixsetup_kappa(kappa, 
                                           expcov=model@sample@covs,
                                           nNode = nNode, 
                                           nGroup = nGroup, 
                                           labels = sampleStats@variables$label,
                                           equal = "kappa" %in% equal, sampletable = sampleStats)
  } else if (type == "cor"){
    # Add rho matrix:
    modMatrices$rho <- matrixsetup_rho(rho, 
                                       expcov=model@sample@covs,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = sampleStats@variables$label,
                                       equal = "rho" %in% equal, sampletable = sampleStats)
    
    if (!corinput){
      # Add SD matrix (ignored if corinput == TRUE):
      modMatrices$SD <- matrixsetup_SD(SD, 
                                       expcov=model@sample@covs,
                                       nNode = nNode, 
                                       nGroup = nGroup, 
                                       labels = sampleStats@variables$label,
                                       equal = "SD" %in% equal, sampletable = sampleStats) 
    }      
  }
  
  
  
  
  
 
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  # if (type == "cov" || type == "prec"){
  #   model@extramatrices <- list(
  #     D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
  #     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
  #     In = as(diag(nNode),"dMatrix") # Identity of dim n
  #   )
  # } else if (type == "chol"){
  #   model@extramatrices <- list(
  #     D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
  #     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
  #     In = as(diag(nNode),"dMatrix"), # Identity of dim n
  #     # C = as(lavaan::lav_matrix_commutation(nNode,nNode),"sparseMatrix")
  #     C = as(lavaan::lav_matrix_commutation(nNode,nNode),"dMatrix")
  #   )
  # } else if (type == "ggm" || type == "cor"){
  #   model@extramatrices <- list(
  #     D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
  #     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
  #     Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
  #     In = as(diag(nNode),"dMatrix"), # Identity of dim n
  #     A = psychonetrics::diagonalizationMatrix(nNode)
  #   )
  # }
  
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
      In = as(diag(nNode),"dMatrix"), # Identity of dim n
      C = as(lavaan::lav_matrix_commutation(nNode,nNode),"dMatrix"),
      Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
      A = psychonetrics::diagonalizationMatrix(nNode)
    )
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Form baseline model:
    if ("omega" %in% equal | "sigma" %in% equal | "kappa" %in% equal){
      equal <- c(equal,"lowertri")
    }
    if (!corinput){
      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "chol",
                                                  lowertri = "diag",
                                                  vars = vars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  meanstructure=meanstructure,
                                                  corinput = corinput,
                                                  ordered = ordered,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats)
    } else {
      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "cor",
                                                  rho = "zero",
                                                  vars = vars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  meanstructure=meanstructure,
                                                  corinput = corinput,
                                                  ordered = ordered,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats)
    }

    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    if (!corinput){
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
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats)
    } else {
      model@baseline_saturated$saturated <- varcov(data,
                                                   type = "cor", 
                                                   lowertri = "full", 
                                                   vars = vars,
                                                   groups = groups,
                                                   covs = covs,
                                                   means = means,
                                                   nobs = nobs,
                                                   missing = missing,
                                                   equal = equal,
                                                   estimator = estimator,
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats)
    }

    
    # if not FIML/FIPML, Treat as computed:
    if (!estimator %in% c("FIML", "FIPML")){
      model@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)      
    }
  }
  
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # Setup PML/FIPML penalization:
  if (estimator %in% c("PML", "FIPML")) {
    model@penalty <- list(lambda = penalty_lambda, alpha = penalty_alpha)
    pen_mats <- if (missing(penalize_matrices)) defaultPenalizeMatrices(model) else penalize_matrices
    model <- penalize(model, matrix = pen_mats, lambda = penalty_lambda, log = FALSE)
    # Baseline/saturated models should use unpenalized estimator:
    base_est <- if (estimator == "FIPML") "FIML" else "ML"
    if (!is.null(model@baseline_saturated$baseline)) {
      model@baseline_saturated$baseline@estimator <- base_est
    }
    if (!is.null(model@baseline_saturated$saturated)) {
      model@baseline_saturated$saturated@estimator <- base_est
    }
  }

  # Return model:
  return(model)
}