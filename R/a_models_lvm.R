# Latent variable model (lvm)
lvm <- function(
  data, # Dataset
  lambda, # only required non-missing matrix
  latent = c("cov","chol","prec","ggm","cor"),
  residual = c("cov","chol","prec","ggm","cor"),

  # Latent matrices:
  sigma_zeta = "full",
  kappa_zeta = "full", # Precision
  omega_zeta = "full", # Partial correlations
  lowertri_zeta = "full", # Cholesky
  delta_zeta = "full", # Used for both ggm and pcor

  # Residual matrices:
  sigma_epsilon = "diag", #
  kappa_epsilon = "diag", # Precision
  omega_epsilon = "zero", # Partial correlations
  lowertri_epsilon = "diag", # Cholesky
  delta_epsilon = "diag", # Used for both ggm and pcor

  # Beta:
  beta = "zero",

  # Mean structure:
  nu,
  nu_eta,
  tau,

  # Identification:
  identify = TRUE,
  identification = c("loadings","variance"),

  # Rest:
  vars, # character indicating the variables Extracted if missing from data - group variable
  ordered = character(0), # character indicating the variables that are ordinal
  latents, # Name of latent varianles
  groups, # deprecated, use groupvar instead
  groupvar, # grouping variable (character column name or vector of group names)
  covs, # alternative covs (array nvar * nvar * ngroup)
  cors, # alternative cors (treated as covs with warning in lvm)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "auto",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  likelihood = c("normal","wishart"), # Gaussian likelihood scaling (ML only). "normal" (default) uses the n denominator; "wishart" uses the n-1 denominator (unbiased S, chisq = (N-1) Fhat, SEs scaled by sqrt(N/(N-1))), matching lavaan likelihood = "wishart".
  fixed_x = character(0), # Exogenous latent variable(s) whose means and mutual (co)variances are fixed to their (free) ML values and excluded from npar/df, matching lavaan fixed.x = TRUE for observed exogenous covariates modelled as single-indicator latents.
  optimizer,
  storedata = FALSE,
  WLS.W,
  covtype = c("choose","ML","UB"),
  corinput, # defaults to FALSE for continuous, TRUE for ordinal
  standardize = c("none","z","quantile"),
  sampleStats,
  verbose=FALSE,
  simplelambdastart = FALSE,
  start = "default",  # "default" (OLS-based), "simple" (old 0.001), or a psychonetrics object
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices,  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
  # Correlation parameterization (latent/residual = "cor") matrices
  # (placed at the end of the signature for backward compatibility):
  rho_zeta = "full",   # Latent correlations
  SD_zeta = "full",    # Latent standard deviations (diagonal)
  rho_epsilon = "zero", # Residual correlations
  SD_epsilon = "diag"  # Residual standard deviations (diagonal)
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
    family = "lvm", caller = "lvm()", estimator = estimator
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

  # Check start:
  start_mod <- NULL
  if (is.character(start)){
    start <- start[1]
    if (!start %in% c("default","simple")){
      stop("'start' can only be 'default', 'simple', or a psychonetrics object")
    }
  } else if (is(start,"psychonetrics")){
    start_mod <- start
    start <- "psychonetrics"
  } else {
    stop("'start' can only be 'default', 'simple', or a psychonetrics object")
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
  # MLR's Huber-White SEs and Yuan-Bentler test need casewise scores, so the
  # raw data must be stored. Force storedata = TRUE for MLR (the MLM family only
  # needs Gamma, which is stored in the sample stats):
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

  # fixed.x exogenous covariates (latent names). Supported for complete-data ML
  # only; the named latents must be exogenous single-indicator latents (the
  # standard way to enter observed exogenous covariates in a psychonetrics SEM).
  fixed_x <- as.character(fixed_x)
  if (length(fixed_x) > 0){
    if (length(ordered) > 0){
      stop("fixed_x is not supported for ordinal data.")
    }
    if (estimator %in% c("FIML","PFIML")){
      stop("fixed_x is not supported for the 'FIML' estimator. Use complete data with estimator = 'ML'.")
    }
  }

  # Gaussian likelihood scaling ("normal" / "wishart"). Wishart is a
  # complete-data ML feature (unbiased S, (N-1) chisq/SE scaling); error
  # informatively for FIML (missing data), ordinal/least-squares and raw
  # time-series input (lavaan likewise restricts likelihood = "wishart" to ML):
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
    experimentalWarning("ordinal data in lvm()")
  }
  if (likelihood == "wishart" && verbose) {
    experimentalWarning("wishart likelihood")
  }

  # Check WLS for ordinal:
  if (length(ordered) > 0 & !estimator %in% c("WLS","DWLS","ULS")){
    stop("Ordinal data is only supported for WLS, DWLS and ULS estimators.")
  }

  # Type:
  latent <- match.arg(latent)
  residual <- match.arg(residual)

  # Identification:
  identification <- match.arg(identification)

  # For ordinal data, force variance identification:
  if (length(ordered) > 0){
    identification <- "variance"
  }

  # WLS weights:
  if (missing(WLS.W)){
    WLS.W <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                    switch(estimator,
                           "WLS" = "full",
                           "ULS" = "identity",
                           "DWLS" = "diag"
                    ))
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

    if (length(ordered) > 0){
      # Ordinal data: pass corinput and meanstructure explicitly
      sampleStats <- samplestats(data = data,
                                 vars = vars,
                                 ordered = ordered,
                                 groups = groups,
                                 covs = covs,
                                 means = means,
                                 nobs = nobs,
                                 missing = ifelse(estimator %in% c("FIML", "PFIML"),"pairwise",missing),
                                 rawts = rawts,
                                 fimldata = estimator %in% c("FIML", "PFIML"),
                                 storedata = storedata,
                                 covtype = covtype,
                                 weightsmatrix = WLS.W,
                                 corinput = TRUE,
                                 meanstructure = FALSE,
                                 verbose = verbose,
                                 standardize = standardize,
                                 bootstrap = bootstrap,
                                 boot_sub = boot_sub,
                                 boot_resample = boot_resample)
    } else {
      sampleStats <- samplestats(data = data,
                                 vars = vars,
                                 groups = groups,
                                 covs = covs,
                                 means = means,
                                 nobs = nobs,
                                 missing = ifelse(estimator %in% c("FIML", "PFIML"),"pairwise",missing),
                                 rawts = rawts,
                                 fimldata = estimator %in% c("FIML", "PFIML"),
                                 storedata = storedata,
                                 covtype = covtype,
                                 weightsmatrix = WLS.W,
                                 verbose = verbose,
                                 standardize = standardize,
                                 bootstrap = bootstrap,
                                 boot_sub = boot_sub,
                                 boot_resample = boot_resample,
                                 likelihood = likelihood)
    }
  }


  # Overwrite corinput:
  corinput <- sampleStats@corinput
  # For continuous data, don't use corinput even if auto-detected (e.g. cor(dat) as covs):
  if (length(ordered) == 0) {
    corinput <- FALSE
    sampleStats@corinput <- FALSE
  }

  # Set meanstructure:
  if (length(ordered) > 0){
    meanstructure <- FALSE
  } else if (corinput){
    # Correlation input: no mean structure
    meanstructure <- FALSE
  } else {
    meanstructure <- TRUE
  }

  # Check some things:
  nNode <- nrow(sampleStats@variables)

  # Generate model object:
  model <- generate_psychonetrics(model = "lvm",sample = sampleStats,computed = FALSE,
                                  equal = equal,identification=identification,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(latent = latent, residual = residual, likelihood = likelihood, fixed_x = fixed_x),
                                  meanstructure = meanstructure,
                                  verbose=verbose)

  # Submodel:
  latentCov <- latent %in% c("cov","chol","cor")
  residCov <- residual %in% c("cov", "chol", "cor")

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

  # Number of thresholds:
  if (length(ordered) > 0){
    nThresh <- sum(sapply(model@sample@thresholds,function(x)sum(sapply(x,length))))
  } else {
    nThresh <- 0
  }

  # Number of means:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  # Add number of observations:
  model@sample@nobs <-
    nNode * (nNode-1) / 2 * nGroup + # Off-diagonal covariances per group
    (!corinput) * nNode * nGroup + # Variances (ignored if correlation matrix is input)
    meanstructure * nMeans + # Means per group
    nThresh # Thresholds
  
  # Stop if lambda is missing:
  if (missing(lambda)){
    stop("'lambda' may not be missing")
  }
  # Handle lambda = "full":
  if (is.character(lambda)){
    if (lambda == "full") {
      if (missing(latents)) {
        stop("'latents' must be specified when lambda = 'full'")
      }
      lambda <- matrix(1, nrow = nNode, ncol = length(latents))
    } else {
      stop("'lambda' must be a numeric matrix, or 'full'")
    }
  }
  
  # Number of latents:
  nLatent <- ncol(lambda)

  # Warn (do not error: may be legitimate with additional constraints) if there
  # are more latent than observed variables; in that case the automatic factor-
  # analytic start values cannot be computed and the model silently falls back
  # to "simple" starting values:
  if (nLatent > nNode){
    warning(paste0("Model has more latent variables (", nLatent, ") than observed variables (", nNode, "). This model may not be identified, and simple starting values will be used."))
  }


  # If latents is not provided, make it:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }
  
  # Model matrices:
  modMatrices <- list()

  # Fix nu (only if mean structure is modeled):
  if (meanstructure){
    modMatrices$nu <- matrixsetup_mu(nu,nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "nu" %in% equal,
                         expmeans = model@sample@means, sampletable = sampleStats, name = "nu")

    # Fix nu_eta
    modMatrices$nu_eta <- matrixsetup_mu(nu_eta,nNode = nLatent,nGroup = nGroup,labels = latents,equal = "nu_eta" %in% equal,
                                      expmeans = lapply(seq_len(nGroup),function(x)rep(0,nLatent)), sampletable = sampleStats, name = "nu_eta")
  }

  # Thresholds:
  if (length(ordered) > 0){
    modMatrices$tau <- matrixsetup_tau(tau, nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,
                                       equal = "tau" %in% equal, sampleThresholds = model@sample@thresholds, sampletable = sampleStats)
  }
  
   # Setup lambda:
  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov=model@sample@covs, nGroup = nGroup, equal = "lambda" %in% equal,
                                           observednames = sampleStats@variables$label, latentnames = latents,
                                           sampletable = sampleStats, identification = identification, simple = simplelambdastart)

  # If nu is entirely fixed but nu_eta has free elements (e.g., latent growth
  # models), all observed means must be absorbed by nu_eta. A start of zero can
  # then be so far off that the optimizer diverges, so instead start nu_eta at
  # the least-squares solution ginv(Lambda) %*% (observed means - fixed nu):
  if (meanstructure){
    nu_anyfree <- any(modMatrices$nu[[1]] != 0)
    nu_eta_free <- modMatrices$nu_eta[[1]] != 0
    if (!nu_anyfree && any(nu_eta_free)){
      for (g in seq_len(nGroup)){
        obsmeans <- as.vector(unlist(model@sample@means[[g]]))
        if (length(obsmeans) == nNode && !any(is.na(obsmeans))){
          lambdastart_g <- as.matrix(modMatrices$lambda$start[,,g])
          nu_eta_ls <- as.vector(MASS::ginv(lambdastart_g) %*% (obsmeans - modMatrices$nu$start[,g]))
          modMatrices$nu_eta$start[nu_eta_free[,g],g] <- nu_eta_ls[nu_eta_free[,g]]
        }
      }
    }
  }

  # Compute the expected latent and residual cov matrices:
  expLatSigma <- lapply(1:nGroup,function(x)matrix(0,nLatent,nLatent))
  expResidSigma <- lapply(1:nGroup,function(x)matrix(0,nNode,nNode))

  # For each group:
  for (g in 1:nGroup){
    expResidSigma[[g]] <- modMatrices$lambda$sigma_epsilon_start[,,g]
    expLatSigma[[g]] <-  modMatrices$lambda$sigma_zeta_start[,,g]
  }

  # Setup beta with OLS-based starting values:
  # beta structure (3D array) is already fixed by fixMatrix inside matrixsetup_beta,
  # but we need it fixed here to determine which elements are non-zero:
  betaFixed <- fixMatrix(beta, nGroup = nGroup, nrows = nLatent, ncols = nLatent,
                         equal = "beta" %in% equal, diag0 = TRUE)

  # Determine if beta has any non-zero structure:
  betaHasStructure <- any(betaFixed != 0)

  if (betaHasStructure && start == "default") {
    # OLS-based beta starting values using sigma_zeta_start (= total Cov(eta) from FA)
    beta_start_list <- vector("list", nGroup)
    use_ols <- TRUE

    # Check overall quality of sigma_zeta_start across all groups:
    for (g in 1:nGroup) {
      S <- as.matrix(expLatSigma[[g]])

      # Guard rail 1: Check if S is near-singular
      ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      if (any(!is.finite(ev)) || min(Re(ev)) / max(Re(ev)) < 1e-8) {
        use_ols <- FALSE
        break
      }
    }

    if (use_ols) {
      for (g in 1:nGroup) {
        S <- as.matrix(expLatSigma[[g]])
        betaMat <- betaFixed[,,g]
        beta_g <- matrix(0, nLatent, nLatent)

        for (i in 1:nLatent) {
          predictors <- which(betaMat[i,] != 0)
          if (length(predictors) == 0) next

          # Extract predictor submatrix and cross-covariance:
          S_pp <- S[predictors, predictors, drop = FALSE]
          s_pi <- S[predictors, i, drop = FALSE]

          # Guard rail 2: Check condition of predictor block
          if (length(predictors) > 1) {
            ev_pp <- eigen(S_pp, symmetric = TRUE, only.values = TRUE)$values
            if (any(!is.finite(ev_pp)) || max(Re(ev_pp)) / min(Re(ev_pp)) > 1e6) {
              use_ols <- FALSE
              break
            }
          }

          # Compute OLS coefficients:
          ols_result <- tryCatch({
            as.numeric(solve(S_pp, s_pi))
          }, error = function(e) NULL)

          if (is.null(ols_result) || any(!is.finite(ols_result))) {
            use_ols <- FALSE
            break
          }

          beta_g[i, predictors] <- ols_result
        }

        if (!use_ols) break
        beta_start_list[[g]] <- beta_g
      }
    }

    # Guard rail 3: Check for extreme coefficients
    if (use_ols) {
      for (g in 1:nGroup) {
        nonzero <- which(betaFixed[,,g] != 0)
        if (any(abs(beta_start_list[[g]][nonzero]) > 5)) {
          use_ols <- FALSE
          break
        }
      }
    }

    # Guard rail 4: Check (I - Beta_start) is well-conditioned
    if (use_ols) {
      for (g in 1:nGroup) {
        IminB <- diag(nLatent) - beta_start_list[[g]]
        IminB_inv <- tryCatch(solve(IminB), error = function(e) NULL)
        if (is.null(IminB_inv) || max(abs(IminB_inv)) > 100) {
          use_ols <- FALSE
          break
        }
      }
    }

    if (use_ols) {
      # Pass OLS starts to matrixsetup_beta:
      modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                            labels = latents, sampletable = sampleStats,
                                            equal = "beta" %in% equal,
                                            start = beta_start_list)

      # Adjust expLatSigma from total Cov(eta) to disturbance Cov(zeta):
      # Sigma_zeta = (I - B) %*% Cov(eta) %*% t(I - B)
      for (g in 1:nGroup) {
        IminB <- diag(nLatent) - beta_start_list[[g]]
        S_total <- as.matrix(expLatSigma[[g]])
        sigma_zeta_adj <- IminB %*% S_total %*% t(IminB)

        # Check if adjusted covariance is positive definite:
        sigma_zeta_adj <- spectralshift(sigma_zeta_adj)
        expLatSigma[[g]] <- sigma_zeta_adj
      }
    } else {
      # Fallback to simple 0.001 starts:
      modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                            labels = latents, sampletable = sampleStats,
                                            equal = "beta" %in% equal)
    }

  } else if (betaHasStructure && start == "psychonetrics") {
    # Extract beta from fitted psychonetrics model:
    beta_est <- getmatrix(start_mod, "beta")
    if (!is.list(beta_est)) {
      beta_est <- list(beta_est)
    }
    modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                          labels = latents, sampletable = sampleStats,
                                          equal = "beta" %in% equal,
                                          start = beta_est)

  } else {
    # Simple starts (0.001) or beta = "zero":
    modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                          labels = latents, sampletable = sampleStats,
                                          equal = "beta" %in% equal)
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
                                           equal = "omega_zeta" %in% equal, sampletable = sampleStats,
                                           beta = modMatrices$beta[[1]])

    # Add delta matrix:
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta, 
                                                name = "delta_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent, 
                                                nGroup = nGroup, 
                                                labels = latents,
                                           equal = "delta_zeta" %in% equal, sampletable = sampleStats,
                                           omegaStart =  modMatrices$omega_zeta$start) 
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
  } else if (latent == "cor"){
    # Add rho matrix:
    modMatrices$rho_zeta <- matrixsetup_rho(rho_zeta,
                                                name = "rho_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                           equal = "rho_zeta" %in% equal, sampletable = sampleStats,
                                           beta = modMatrices$beta[[1]])

    # Add SD matrix:
    modMatrices$SD_zeta <- matrixsetup_SD(SD_zeta,
                                                name = "SD_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                           equal = "SD_zeta" %in% equal, sampletable = sampleStats,
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
                                                equal = "delta_epsilon" %in% equal, sampletable = sampleStats,
                                                omegaStart =  modMatrices$omega_epsilon$start) 
  } else if (residual == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_epsilon <- matrixsetup_kappa(kappa_epsilon,
                                                   name = "kappa_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode,
                                                nGroup = nGroup,
                                                labels = sampleStats@variables$label,
                                                equal = "kappa_epsilon" %in% equal, sampletable = sampleStats)
  } else if (residual == "cor"){
    # Add rho matrix:
    modMatrices$rho_epsilon <- matrixsetup_rho(rho_epsilon,
                                                   name = "rho_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode,
                                                nGroup = nGroup,
                                                labels = sampleStats@variables$label,
                                                equal = "rho_epsilon" %in% equal, sampletable = sampleStats)

    # Add SD matrix:
    modMatrices$SD_epsilon <- matrixsetup_SD(SD_epsilon,
                                                   name = "SD_epsilon",
                                                expcov=expResidSigma,
                                                nNode = nNode,
                                                nGroup = nGroup,
                                                labels = sampleStats@variables$label,
                                                equal = "SD_epsilon" %in% equal, sampletable = sampleStats)
  }
  
  
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  # When corinput=TRUE, fix diagonal elements of sigma_epsilon (they are derived from diag(sigma)=1):
  if (corinput){
    # Find the diagonal parameters of residual matrices:
    diagFixMatrices <- c("sigma_epsilon", "kappa_epsilon", "lowertri_epsilon", "delta_epsilon", "SD_epsilon")
    for (matName in diagFixMatrices){
      diagRows <- pars$partable$matrix == matName & pars$partable$row == pars$partable$col
      if (any(diagRows)){
        pars$partable$fixed[diagRows] <- TRUE
        pars$partable$par[diagRows] <- 0
        # Set starting value: for sigma_epsilon cov, start at 0.5; for others use current est
        # The actual values will be derived in the implied function
      }
    }
  }

  # fixed.x: fix the exogenous latent block (their sigma_zeta (co)variances and
  # the observed indicators' intercepts) to the sample moments and remove their
  # statistics from the df count. The named exogenous variables must be
  # single-indicator latents (the standard way to enter observed exogenous
  # covariates into a psychonetrics SEM): each must be a latent whose lambda
  # column has exactly one non-zero (unit) loading.
  if (length(fixed_x) > 0){
    bad <- setdiff(fixed_x, latents)
    if (length(bad) > 0){
      stop("fixed_x must name exogenous *latent* variables (the single-indicator latents carrying the observed covariates). Not found among the latents: ",
           paste(bad, collapse = ", "), ". The latents are: ", paste(latents, collapse = ", "), ".")
    }
    lat_idx <- match(fixed_x, latents)
    obs_idx <- integer(length(lat_idx))
    for (k in seq_along(lat_idx)){
      col <- lambda[, lat_idx[k]]
      nz <- which(col != 0)
      if (length(nz) != 1){
        stop("fixed_x latent '", fixed_x[k], "' is not a single-indicator latent (its lambda column must have exactly one non-zero loading). fixed_x in lvm() is only supported for observed exogenous covariates entered as single-indicator latents.")
      }
      obs_idx[k] <- nz
    }

    pars$partable <- apply_fixed_x_partable_lvm(
      partable = pars$partable,
      lat_idx = lat_idx,
      obs_idx = obs_idx,
      sample_covs = model@sample@covs,
      sample_means = model@sample@means,
      group_ids = model@sample@groups$id,
      meanstructure = meanstructure
    )

    # Remove the x-block statistics from the df count:
    p_x <- length(fixed_x)
    model@sample@nobs <- model@sample@nobs - fixed_x_nstat_drop(p_x, nGroup, meanstructure)

    # Record the observed indicator indices so the conditional-logl adjustment
    # (fixed_x_marginal_loglik) can find the x-block among the observed variables:
    model@types$fixed_x_obs_idx <- obs_idx

    if (verbose) experimentalWarning("fixed.x exogenous covariates")
  }

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
    In = as(diag(nNode),"dMatrix"), # Identity of dim n
    Inlatent = as(diag(nLatent),"dMatrix"),
    C = commutationMatrix(nNode, nLatent),
    Cbeta = commutationMatrix(nLatent, nLatent),
    C_chol = commutationMatrix(nNode, nNode),
    A = psychonetrics::diagonalizationMatrix(nNode),
    Aeta = psychonetrics::diagonalizationMatrix(nLatent)
  )
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){

    if (!corinput){
      # Form baseline model:
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
                                               meanstructure = meanstructure,
                                               corinput = corinput,
                                               ordered = ordered,
                                               baseline_saturated = FALSE, sampleStats = sampleStats, likelihood = likelihood)
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
                                                  meanstructure = meanstructure,
                                                  corinput = corinput,
                                                  ordered = ordered,
                                                  baseline_saturated = FALSE, sampleStats = sampleStats, likelihood = likelihood)
    }

    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)


    ### Saturated model ###
    # No `equal = equal`: the saturated reference is unconstrained per
    # group; cross-group equality belongs to the target/baseline only.
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
             meanstructure = meanstructure,
             corinput = corinput,
             ordered = ordered,
             baseline_saturated = FALSE, sampleStats = sampleStats, likelihood = likelihood)
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
             meanstructure = meanstructure,
             corinput = corinput,
             ordered = ordered,
             baseline_saturated = FALSE, sampleStats = sampleStats, likelihood = likelihood)
    }

    # if not FIML/PFIML, Treat as computed:
    if (!estimator %in% c("FIML", "PFIML")){
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

  # Store robust ML configuration (MLM/MLMV/MLMVS/MLR) and propagate it to the
  # baseline/saturated models so the scaled baseline statistic uses the same
  # correction (needed for the robust incremental fit indices):
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
