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
  groups, # deprecated, use groupvar instead
  groupvar, # grouping variable (character column name or vector of group names)
  covs, # alternative covs (array nvar * nvar * ngroup)
  cors, # alternative cors (treated as covs with corinput=TRUE for varcov)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "auto",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  likelihood = c("normal","wishart"), # Gaussian likelihood scaling (ML/FIML only). "normal" (default) uses the n denominator; "wishart" uses the n-1 denominator (unbiased S, chisq = (N-1) Fhat, SEs scaled by sqrt(N/(N-1))), matching lavaan likelihood = "wishart".
  fixed_x = character(0), # Exogenous variables whose means and mutual (co)variances are fixed to their sample values (excluded from npar/df), matching lavaan fixed.x = TRUE. Only for type = "cov".
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
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
){
  # Standardize input arguments:
  si <- standardize_input(
    data = if(missing(data)) NULL else data,
    covs = if(missing(covs)) NULL else covs,
    cors = if(missing(cors)) NULL else cors,
    nobs = if(missing(nobs)) NULL else nobs,
    corinput = if(missing(corinput)) NULL else corinput,
    groups = if(missing(groups)) NULL else groups,
    groupvar = if(missing(groupvar)) NULL else groupvar,
    family = "varcov", caller = "varcov()", estimator = estimator
  )
  # Only overwrite when standardize_input actually resolved a value,
  # so that the original missing() state is preserved for downstream code:
  if (!is.null(si$data)) data <- si$data
  if (!is.null(si$covs)) covs <- si$covs
  if (!is.null(si$nobs)) nobs <- si$nobs
  if (!is.null(si$corinput)) corinput <- si$corinput
  if (!is.null(si$groups)) groups <- si$groups

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

  # Robust ML estimators (MLM/MLMV/MLMVS/MLR) map to estimator = "ML" plus a
  # robust SE/test configuration (Phase 1: complete data only):
  .robust_resolved <- resolve_robust_estimator(estimator)
  estimator <- .robust_resolved$estimator
  robust_cfg <- .robust_resolved$robust
  if (isTRUE(nzchar(robust_cfg$se)) && robust_cfg$label == "MLR"){
    storedata <- TRUE
  }

  if (estimator == "default"){
    if (length(ordered) > 0){
      estimator <- "DWLS"
    } else {
      estimator <- "ML"
    }
  }

  # Gaussian likelihood scaling ("normal" / "wishart"). The wishart scaling is a
  # complete-data ML feature: it requires the n-1 (unbiased) sample covariance and
  # the (N-1) chisq/SE scaling, which are only defined for complete-data ML (as in
  # lavaan, which likewise restricts likelihood = "wishart" to ML). Error
  # informatively for FIML / PFIML (missing data) and the least-squares
  # estimators:
  likelihood <- match.arg(likelihood)
  if (likelihood == "wishart"){
    if (estimator %in% c("FIML","PFIML")){
      stop("likelihood = 'wishart' is not supported for the 'FIML' estimator (it requires complete-data maximum likelihood). Use likelihood = 'normal' (the default).")
    }
    if (estimator %in% c("WLS","DWLS","ULS","PML")){
      stop("likelihood = 'wishart' is only supported for maximum-likelihood (estimator = 'ML') estimation, not for ", estimator, ".")
    }
    if (length(ordered) > 0){
      stop("likelihood = 'wishart' is not supported for ordinal data.")
    }
    if (isTRUE(rawts)){
      stop("likelihood = 'wishart' is not supported for raw time-series input.")
    }
  }

  # Experimental warnings:
  if (length(ordered) > 0) {
    experimentalWarning("ordinal data in varcov()")
  }
  if (likelihood == "wishart" && verbose) {
    experimentalWarning("wishart likelihood")
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

  # fixed.x exogenous covariates. Only supported for the covariance
  # parameterisation (type = "cov"), where the x-x block lives directly in the
  # 'sigma' matrix; the other parameterisations (ggm/prec/chol/cor) do not
  # expose the x-x covariance block as free parameters.
  fixed_x <- as.character(fixed_x)
  if (length(fixed_x) > 0){
    if (type != "cov"){
      stop("fixed_x is currently only supported for varcov(type = 'cov').")
    }
    if (length(ordered) > 0){
      stop("fixed_x is not supported for ordinal data.")
    }
    if (estimator %in% c("FIML","PFIML")){
      stop("fixed_x is not supported for the 'FIML' estimator. Use complete data with estimator = 'ML'.")
    }
  }

  # Set meanstructure:
  if (missing(meanstructure)){
    meanstructure <- (!missing(data) || !missing(means))
  }

  # Check FIML:
  if (!missing(data) && !meanstructure && estimator %in% c("FIML", "PFIML")){
    stop("meanstructure = FALSE is not yet supported for 'FIML'/'PFIML' estimator")
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
        if (likelihood == "wishart") {
          stop("likelihood = 'wishart' is not supported with missing data (it requires complete-data maximum likelihood). Use listwise deletion or likelihood = 'normal'.")
        }
        if (estimator == "ML") {
          estimator <- "FIML"
        } else if (estimator == "PML") {
          estimator <- "PFIML"
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
                               missing = ifelse(estimator %in% c("FIML", "PFIML"),"pairwise",missing),
                               rawts = rawts,
                               fimldata = estimator %in% c("FIML", "PFIML"),
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
                               boot_resample = boot_resample,
                               likelihood = likelihood)
  }

  # Overwrite corinput:
  corinput <- sampleStats@corinput

  # Meanstructure is not supported in combination with correlation matrix input
  # for maximum likelihood estimation (corinput = TRUE implies standardized data,
  # so a meanstructure is not meaningful):
  if (corinput && meanstructure && estimator %in% c("ML","PML")){
    stop("meanstructure = TRUE is not supported in combination with corinput = TRUE for (penalized) maximum likelihood estimation. Use meanstructure = FALSE, as correlation matrix input implies standardized data.")
  }

  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Generate model object:
  # Validate fixed_x against the variable names (now that sample stats exist):
  if (length(fixed_x) > 0){
    varlabels <- sampleStats@variables$label
    bad <- setdiff(fixed_x, varlabels)
    if (length(bad) > 0){
      stop("fixed_x variable(s) not found among the model variables: ", paste(bad, collapse = ", "))
    }
    if (length(fixed_x) >= length(varlabels)){
      stop("fixed_x cannot include all variables: at least one endogenous variable must remain.")
    }
    if (verbose) experimentalWarning("fixed.x exogenous covariates")
  }

  model <- generate_psychonetrics(model = "varcov",sample = sampleStats,computed = FALSE,
                                  equal = equal,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(y = type, likelihood = likelihood, fixed_x = fixed_x),
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

  # fixed.x: the pure x-block (x means + x-x (co)variances) is conditioned on
  # rather than modelled, so its statistics leave the df count:
  if (length(fixed_x) > 0){
    p_x <- length(fixed_x)
    model@sample@nobs <- model@sample@nobs - fixed_x_nstat_drop(p_x, nGroup, meanstructure)
  }

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


  # fixed.x: fix the x-block of the parameter table (x means in 'mu' and the x-x
  # covariance block in 'sigma') to the sample moments and renumber the free
  # parameters. The cross (x-y) covariances and the y-y block remain free, so
  # the model conditions on x while every endogenous moment stays estimated:
  if (length(fixed_x) > 0){
    x_idx <- match(fixed_x, sampleStats@variables$label)
    pars$partable <- apply_fixed_x_partable(
      partable = pars$partable,
      x_idx = x_idx,
      cov_matrix = "sigma",
      mean_matrix = if (meanstructure) "mu" else NA_character_,
      sample_covs = model@sample@covs,
      sample_means = model@sample@means,
      group_ids = model@sample@groups$id
    )
  }

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
      C = commutationMatrix(nNode,nNode),
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
                                                  baseline_saturated = FALSE,sampleStats=sampleStats, likelihood = likelihood)
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
                                                  baseline_saturated = FALSE,sampleStats=sampleStats, likelihood = likelihood)
    }

    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    # The saturated reference model is intentionally left unconstrained
    # (no `equal = equal`): the LR-based fit measures must reference a
    # model in which each group's covariance/threshold structure is fully
    # free. Cross-group equality is a property of the target/baseline,
    # not the saturated. (Was a long-standing bug pre-0.15.10.)
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
                                                   estimator = estimator,
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats, likelihood = likelihood)
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
                                                   estimator = estimator,
                                                   meanstructure=meanstructure,
                                                   corinput = corinput,
                                                   ordered = ordered,
                                                   baseline_saturated = FALSE,sampleStats=sampleStats, likelihood = likelihood)
    }

    
    # if not FIML/PFIML, Treat as computed:
    if (!estimator %in% c("FIML", "PFIML")){
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

  # Setup PML/PFIML penalization:
  if (estimator %in% c("PML", "PFIML")) {
    model@penalty <- list(lambda = penalty_lambda, alpha = penalty_alpha)
    pen_mats <- if (missing(penalize_matrices)) defaultPenalizeMatrices(model) else penalize_matrices
    model <- penalize(model, matrix = pen_mats, lambda = penalty_lambda, log = FALSE)
    # Baseline/saturated models should use unpenalized estimator:
    base_est <- if (estimator == "PFIML") "FIML" else "ML"
    if (!is.null(model@baseline_saturated$baseline)) {
      model@baseline_saturated$baseline@estimator <- base_est
    }
    if (!is.null(model@baseline_saturated$saturated)) {
      model@baseline_saturated$saturated@estimator <- base_est
    }
  }

  # Store robust ML configuration (MLM/MLMV/MLMVS/MLR) and propagate to the
  # baseline/saturated models (for the robust incremental fit indices):
  if (length(robust_cfg) > 0 && .hasSlot(model, "robust")){
    model@robust <- robust_cfg
    if (!is.null(model@baseline_saturated$baseline) &&
        is(model@baseline_saturated$baseline, "psychonetrics") &&
        .hasSlot(model@baseline_saturated$baseline, "robust")){
      model@baseline_saturated$baseline@robust <- robust_cfg
    }
    if (!is.null(model@baseline_saturated$saturated) &&
        is(model@baseline_saturated$saturated, "psychonetrics") &&
        .hasSlot(model@baseline_saturated$saturated, "robust")){
      model@baseline_saturated$saturated@robust <- robust_cfg
    }
  }

  # Return model:
  return(model)
}
