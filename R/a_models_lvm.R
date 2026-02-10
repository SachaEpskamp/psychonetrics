# Latent variable model (lvm)
lvm <- function(
  data, # Dataset
  lambda, # only required non-missing matrix
  latent = c("cov","chol","prec","ggm"), # Maybe add cor at some point, but not now
  residual = c("cov","chol","prec","ggm"),

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
  covtype = c("choose","ML","UB"),
  standardize = c("none","z","quantile"),
  sampleStats,
  corinput,
  verbose=FALSE,
  simplelambdastart = FALSE,
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
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
    } else if (!missing(corinput) && corinput){
      estimator <- "WLS"
    } else {
      estimator <- "ML"
    }
  }

  # Experimental warnings:
  if (length(ordered) > 0) {
    experimentalWarning("ordinal data in lvm()")
  }

  # Check WLS for ordinal:
  if (length(ordered) > 0 & !estimator %in% c("WLS","DWLS","ULS")){
    stop("Ordinal data is only supported for WLS, DWLS and ULS estimators.")
  }

  # Check WLS for corinput:
  if (!missing(corinput) && corinput && !estimator %in% c("WLS","DWLS","ULS")){
    stop("corinput = TRUE is only supported for WLS, DWLS and ULS estimators.")
  }

  # Type:
  latent <- match.arg(latent)
  residual <- match.arg(residual)

  # Identification:
  identification <- match.arg(identification)

  # For ordinal data, force variance identification and corinput:
  if (length(ordered) > 0){
    identification <- "variance"
    if (!missing(corinput) && !corinput){
      stop("corinput must be TRUE for ordinal data.")
    }
  }

  # For continuous corinput, force variance identification:
  if (!missing(corinput) && corinput && length(ordered) == 0){
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
      # Continuous data: pass corinput if specified
      if (!missing(corinput) && corinput){
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
                                   boot_resample = boot_resample)
      }
    }
  }


  # Overwrite corinput:
  corinput <- sampleStats@corinput

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
                                  rawts = rawts, types = list(latent = latent, residual = residual),
                                  meanstructure = meanstructure,
                                  verbose=verbose)

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
  
  # Setup beta:
  modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup, labels = latents, sampletable = sampleStats, equal = "beta" %in% equal)
  
  # Compute the expected latent and residual cov matrices:
  expLatSigma <- lapply(1:nGroup,function(x)matrix(0,nLatent,nLatent))
  expResidSigma <- lapply(1:nGroup,function(x)matrix(0,nNode,nNode))
  
  # For each group:
  for (g in 1:nGroup){
    expResidSigma[[g]] <- modMatrices$lambda$sigma_epsilon_start[,,g]
    expLatSigma[[g]] <-  modMatrices$lambda$sigma_zeta_start[,,g]
    
      
    # # Current cov estimate:
    # curcov <- as.matrix(sampleStats@covs[[g]])
    # 
    # 
    # # Cur loadings:
    # curLambda <- modMatrices$lambda$start[,,g]
    # if (!is.matrix(curLambda)){
    #   curLambda <- as.matrix(curLambda)
    # }
    # 
    # ### Run factor start ###:
    # # facStart <- factorstart(curcov, modMatrices$lambda[[1]][,,1])
    # 
    # 
    # ###
    # 
    # 
    # # Residual variances, let's start by putting the vars on 1/4 times the observed variances:
    # Theta <- diag(diag(curcov)/4)
    # # Theta <- modMatrices$lambda$thetaStart[,,g]
    # fa <- factanal(covmat = curcov, factors = nLatent, covar = TRUE)
    # 
    # fa <- fa(r = curcov, nfactors = nLatent, rotate = "promax", covar = TRUE)
    # 
    # Theta <- diag(fa$uniquenesses)
    # 
    # # Check if this is positive definite:
    # ev <- eigen(curcov - Theta)$values
    # 
    # # Shrink until it is positive definite:
    # loop <- 0
    # repeat{
    #   ev <- eigen(curcov - Theta)$values
    #   if (loop == 100){
    #     # give up...
    #     
    #     Theta <- diag(nrow(Theta))
    #     break
    #   }
    #   if (all(ev>0)){
    #     break
    #   }
    #   Theta <- Theta * 0.9
    #   loop <- loop + 1
    # }
    # 
    # # Expected residual sigma:
    # expResidSigma[[g]] <- Theta
    # 
    # # This means that the factor-part is expected to be:
    # factorPart <- curcov - Theta
    # 
    # # Let's take a pseudoinverse:
    # inv <- corpcor::pseudoinverse(kronecker(curLambda,curLambda))
    # 
    # # And obtain psi estimate:
    # expLatSigma[[g]] <- matrix(inv %*% as.vector(factorPart),nLatent,nLatent)
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
  }
  
  
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  # When corinput=TRUE, fix diagonal elements of sigma_epsilon (they are derived from diag(sigma)=1):
  if (corinput){
    # Find the diagonal parameters of residual matrices:
    diagFixMatrices <- c("sigma_epsilon", "kappa_epsilon", "lowertri_epsilon", "delta_epsilon")
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
    C = as(lavaan::lav_matrix_commutation(nNode, nLatent),"dMatrix"),
    Cbeta = as(lavaan::lav_matrix_commutation(nLatent, nLatent),"dMatrix"),
    C_chol = as(lavaan::lav_matrix_commutation(nNode, nNode),"dMatrix"),
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
                                               baseline_saturated = FALSE, sampleStats = sampleStats)
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
                                                  baseline_saturated = FALSE, sampleStats = sampleStats)
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
             meanstructure = meanstructure,
             corinput = corinput,
             ordered = ordered,
             baseline_saturated = FALSE, sampleStats = sampleStats)
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
             meanstructure = meanstructure,
             corinput = corinput,
             ordered = ordered,
             baseline_saturated = FALSE, sampleStats = sampleStats)
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

  # Return model:
  return(model)
}