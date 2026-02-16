# Meta-analytic VAR(1) model (one-stage random effects)
meta_var1 <- function(
  data,           # Raw time-series data (use with studyvar)
  covs,           # Pre-computed list of Toeplitz covariance matrices (alternative to data)
  nobs,           # Sample sizes per study (required with covs)

  # Study/time structure:
  vars,           # Character vector of variable names
  studyvar,       # Column name identifying studies
  idvar,          # Subject ID within study (for collating time series)
  dayvar,         # Day variable (passed to tsData)
  beepvar,        # Beep variable (passed to tsData)

  # VAR(1) structure:
  contemporaneous = c("cov","chol","prec","ggm"),
  beta = "full",

  # Contemporaneous effects:
  omega_zeta = "full",
  delta_zeta = "full",
  kappa_zeta = "full",
  sigma_zeta = "full",
  lowertri_zeta = "full",

  # Random effects setup:
  randomEffects = c("chol","cov","prec","ggm","cor"),
  sigma_randomEffects = "full",
  kappa_randomEffects = "full",
  omega_randomEffects = "full",
  lowertri_randomEffects = "full",
  delta_randomEffects = "full",
  rho_randomEffects = "full",
  SD_randomEffects = "full",

  # V matrix estimation:
  Vmats,
  Vmethod = c("individual","pooled"),
  Vestimation = c("averaged","per_study"),

  # Some extra stuff:
  baseline_saturated = TRUE,
  optimizer,
  estimator = c("FIML","ML"),

  sampleStats,
  verbose = FALSE,
  bootstrap = FALSE,
  boot_sub,
  boot_resample
){

  message(paste0("Note: 'meta_var1()' is experimental in psychonetrics ",
                 utils::packageVersion("psychonetrics"),
                 ". Please report any unexpected behavior to https://github.com/SachaEpskamp/psychonetrics/issues"))

  contemporaneous <- match.arg(contemporaneous)
  randomEffects <- match.arg(randomEffects)
  Vmethod <- match.arg(Vmethod)
  Vestimation <- match.arg(Vestimation)
  estimator <- match.arg(estimator)

  # Handle idvar/studyvar logic:
  if (!missing(data) && !is.null(data)){
    # Data mode:
    if (missing(studyvar) || is.null(studyvar)){
      if (!missing(idvar) && !is.null(idvar)){
        warning("'studyvar' not supplied; using 'idvar' as study variable. Set 'studyvar' explicitly to suppress this warning.")
        studyvar <- idvar
        idvar <- NULL
      } else {
        stop("'studyvar' is required when using raw data input.")
      }
    }
  }

  # Resolve missing optionals:
  if (missing(idvar)) idvar <- NULL
  if (missing(dayvar)) dayvar <- NULL
  if (missing(beepvar)) beepvar <- NULL
  if (missing(vars)) vars <- NULL

  #### Process input data ####
  if (!missing(data) && !is.null(data)){
    data <- as.data.frame(data)

    # Check studyvar exists:
    if (!studyvar %in% names(data)){
      stop(paste0("'studyvar' column '", studyvar, "' not found in data."))
    }

    # Get study IDs:
    studies <- unique(data[[studyvar]])
    nStudy <- length(studies)

    # If vars not specified, infer from data (exclude structural vars):
    structVars <- c(studyvar, idvar, dayvar, beepvar)
    if (is.null(vars)){
      vars <- names(data)[!names(data) %in% structVars]
    }
    nNode <- length(vars)

    # For each study, compute Toeplitz covariance and extract blocks:
    studyCovs <- list()
    sampleSizes <- numeric(nStudy)

    for (i in seq_len(nStudy)){
      studyData <- data[data[[studyvar]] == studies[i], , drop = FALSE]
      # Remove studyvar column:
      studyData <- studyData[, names(studyData) != studyvar, drop = FALSE]

      # Run tsData to create lagged augmented data:
      augData <- tsData(studyData, vars = vars, beepvar = beepvar, dayvar = dayvar, idvar = idvar)

      # Compute covariance (use pairwise complete):
      covMat <- cov(augData, use = "pairwise.complete.obs")

      # Store:
      studyCovs[[i]] <- covMat
      sampleSizes[i] <- nrow(augData)
    }

    covs <- studyCovs
    nobs <- sampleSizes

  } else if (!missing(covs) && !is.null(covs)){
    # Pre-computed covariance matrices mode:
    if (missing(nobs) || is.null(nobs)){
      stop("'nobs' is required when using pre-computed covariance matrices.")
    }

    nStudy <- length(covs)
    if (length(nobs) == 1) nobs <- rep(nobs, nStudy)
    if (length(nobs) != nStudy){
      stop("Length of 'nobs' does not match number of covariance matrices.")
    }

    # Get nNode from first matrix:
    nVarFull <- ncol(covs[[1]])
    if (nVarFull %% 2 != 0){
      stop("Covariance matrices must be Toeplitz matrices with an even number of variables.")
    }
    nNode <- nVarFull / 2

    # Get var names:
    if (is.null(vars)){
      if (!is.null(colnames(covs[[1]]))){
        # Use names from the endogenous (second) block:
        allNames <- colnames(covs[[1]])
        vars <- allNames[(nNode + 1):(2 * nNode)]
        # Remove _lag suffix if present:
        vars <- sub("_lag[0-9]+$", "", vars)
      } else {
        vars <- paste0("V", seq_len(nNode))
      }
    }

    sampleSizes <- nobs

  } else {
    stop("Either 'data' (with 'studyvar') or 'covs' (with 'nobs') must be supplied.")
  }

  # Extract Sigma0 and Sigma1 blocks from each study's Toeplitz matrix:
  Sigma0_list <- list()
  Sigma1_list <- list()

  for (i in seq_along(covs)){
    fullCov <- as.matrix(covs[[i]])
    Sigma0_list[[i]] <- fullCov[nNode + (1:nNode), nNode + (1:nNode)]
    Sigma1_list[[i]] <- fullCov[nNode + (1:nNode), 1:nNode]
  }

  # Build meta-level data: each study is a row with [vech(Sigma0), vec(Sigma1)]
  # Labels for the meta-level variables:
  sigma0Labels <- paste0(vars[row(matrix(0, nNode, nNode))[lower.tri(matrix(0, nNode, nNode), diag=TRUE)]],
                         " ~~ ",
                         vars[col(matrix(0, nNode, nNode))[lower.tri(matrix(0, nNode, nNode), diag=TRUE)]])
  sigma1Labels <- paste0(vars[row(matrix(0, nNode, nNode))],
                         " -> ",
                         vars[col(matrix(0, nNode, nNode))])

  covvars <- c(sigma0Labels, sigma1Labels)
  nCov <- length(covvars)

  metaData <- do.call(rbind, lapply(seq_along(covs), function(i){
    c(Vech(Sigma0_list[[i]], diag = TRUE), as.vector(Sigma1_list[[i]]))
  }))
  metaData <- as.data.frame(metaData)
  names(metaData) <- covvars

  #### Compute V matrices (sampling error approximation) ####
  if (!missing(Vmats)){
    avgVmat <- Reduce("+", Vmats) / length(Vmats)
  } else {
    if (verbose){
      message("Computing sampling error approximation...")
    }

    if (Vmethod == "individual"){

      # Selection matrix: maps vech(full Toeplitz) -> [vech(Sigma0), vec(Sigma1)]
      # We need to figure out which elements of vech(full) correspond to our elements

      Vmats <- lapply(seq_along(covs), function(i){
        fullCov <- as.matrix(covs[[i]])
        nFull <- nrow(fullCov)
        k <- solve_symmetric_cpp_matrixonly(fullCov)
        Dfull <- duplicationMatrix(nFull)
        v <- 0.5 * sampleSizes[i] * t(Dfull) %*% (k %x% k) %*% Dfull
        vcov_full <- solve_symmetric_cpp_matrixonly(as.matrix(v))

        # Now we need to select the elements corresponding to our [vech(Sigma0), vec(Sigma1)]
        # Build index mapping:
        # vech(full) indices for Sigma0 block: lower.tri of block [nNode+(1:nNode), nNode+(1:nNode)]
        # vec indices for Sigma1 block: all elements of block [nNode+(1:nNode), 1:nNode]
        dumFull <- matrix(0, nFull, nFull)
        dumFull[lower.tri(dumFull, diag = TRUE)] <- seq_len(nFull * (nFull + 1) / 2)

        # Indices in vech(full) for vech(Sigma0):
        sigma0Block <- dumFull[nNode + (1:nNode), nNode + (1:nNode)]
        sigma0Inds <- sigma0Block[lower.tri(sigma0Block, diag = TRUE)]

        # Indices in vech(full) for vec(Sigma1):
        # Sigma1 = full[nNode+(1:nNode), 1:nNode]
        # These are below-diagonal in the full matrix, so they're in the lower triangle
        sigma1Block <- dumFull[nNode + (1:nNode), 1:nNode]
        sigma1Inds <- as.vector(sigma1Block)

        # Selection matrix:
        allInds <- c(sigma0Inds, sigma1Inds)
        nSel <- length(allInds)
        Sel <- sparseMatrix(i = seq_len(nSel), j = allInds,
                            dims = c(nSel, nFull * (nFull + 1) / 2))
        Sel <- as(Sel, "dMatrix")

        # Extract relevant sub-block of vcov:
        res <- as.matrix(Sel %*% vcov_full %*% t(Sel))
        return(0.5 * (res + t(res)))
      })

      avgVmat <- Reduce("+", Vmats) / Reduce("+", lapply(Vmats, function(x) x != 0))

    } else if (Vmethod == "pooled"){
      # Pool all Toeplitz matrices and compute V from pooled estimates:
      # Fit a single varcov model with equal covariances across studies:
      pooledCov <- Reduce("+", lapply(seq_along(covs), function(i) sampleSizes[i] * covs[[i]])) /
        sum(sampleSizes)
      nFull <- nrow(pooledCov)
      k <- solve_symmetric_cpp_matrixonly(as.matrix(pooledCov))
      Dfull <- duplicationMatrix(nFull)
      v <- 0.5 * mean(sampleSizes) * t(Dfull) %*% (k %x% k) %*% Dfull
      vcov_full <- solve_symmetric_cpp_matrixonly(as.matrix(v))

      # Build selection matrix (same as individual):
      dumFull <- matrix(0, nFull, nFull)
      dumFull[lower.tri(dumFull, diag = TRUE)] <- seq_len(nFull * (nFull + 1) / 2)

      sigma0Block <- dumFull[nNode + (1:nNode), nNode + (1:nNode)]
      sigma0Inds <- sigma0Block[lower.tri(sigma0Block, diag = TRUE)]
      sigma1Block <- dumFull[nNode + (1:nNode), 1:nNode]
      sigma1Inds <- as.vector(sigma1Block)
      allInds <- c(sigma0Inds, sigma1Inds)

      nSel <- length(allInds)
      Sel <- sparseMatrix(i = seq_len(nSel), j = allInds,
                          dims = c(nSel, nFull * (nFull + 1) / 2))
      Sel <- as(Sel, "dMatrix")

      avgVmat <- as.matrix(Sel %*% vcov_full %*% t(Sel))
      avgVmat <- 0.5 * (avgVmat + t(avgVmat))

      # Per-study V matrices scaled by sample size:
      Vmats <- lapply(sampleSizes, function(n) mean(sampleSizes) / n * avgVmat)
    }
  }

  #### Obtain sample stats ####
  if (missing(sampleStats)){
    sampleStats <- samplestats(data = metaData,
                               vars = covvars,
                               missing = "pairwise",
                               fimldata = estimator == "FIML",
                               storedata = FALSE,
                               meanstructure = TRUE,
                               verbose = verbose,
                               fullFIML = (Vestimation == "per_study"),
                               bootstrap = bootstrap,
                               boot_sub = boot_sub,
                               boot_resample = boot_resample)
  }

  # The meta-level data is not corinput:
  sampleStats@corinput <- FALSE

  #### Generate model object ####
  model <- generate_psychonetrics(model = "meta_var1",
                                  submodel = switch(contemporaneous,
                                                    "prec" = "meta_gvar",
                                                    "ggm" = "meta_gvar",
                                                    "chol" = "meta_var",
                                                    "cov" = "meta_var"),
                                  sample = sampleStats,
                                  computed = FALSE,
                                  optimizer = defaultoptimizer(),
                                  estimator = estimator,
                                  distribution = "Gaussian",
                                  types = list(zeta = contemporaneous, randomEffects = randomEffects),
                                  meanstructure = TRUE,
                                  verbose = verbose)

  # Number of groups:
  nGroup <- 1  # Multi-group not yet supported

  # Add number of observations:
  nMeans <- sum(sapply(model@sample@means, function(x) sum(!is.na(x))))
  model@sample@nobs <-
    nCov * (nCov - 1) / 2 * nGroup +  # Covariances
    nCov * nGroup +                     # Variances
    nMeans                              # Means

  #### VAR(1) Model matrices ####
  modMatrices <- list()

  # Compute starting values from pooled covariances:
  # Use mean across studies:
  meanSigma0 <- Reduce("+", Sigma0_list) / nStudy
  meanSigma1 <- Reduce("+", Sigma1_list) / nStudy

  # Prior estimate for beta (S1 * S0^{-1}):
  S0inv <- solve_symmetric_cpp_matrixonly(as.matrix(spectralshift(meanSigma0)))
  betaEst <- list(as.matrix(meanSigma1 %*% S0inv))

  modMatrices$beta <- matrixsetup_beta(beta,
                                        name = "beta",
                                        nNode = nNode,
                                        nGroup = nGroup,
                                        labels = vars,
                                        equal = FALSE,
                                        sampletable = sampleStats,
                                        start = betaEst,
                                        onlyStartSign = FALSE)

  # Prior guess for contemporaneous covariances (Schur complement):
  contCovEst <- list(as.matrix(spectralshift(meanSigma0 - meanSigma1 %*% S0inv %*% t(meanSigma1))))

  # Fill in contemporaneous structure:
  if (contemporaneous == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta, name = "sigma_zeta",
                                                 expcov = contCovEst,
                                                 nNode = nNode,
                                                 nGroup = nGroup,
                                                 labels = vars,
                                                 equal = FALSE, sampletable = sampleStats)
  } else if (contemporaneous == "chol"){
    modMatrices$lowertri_zeta <- matrixsetup_lowertri(lowertri_zeta, name = "lowertri_zeta",
                                                       expcov = contCovEst,
                                                       nNode = nNode,
                                                       nGroup = nGroup,
                                                       labels = vars,
                                                       equal = FALSE, sampletable = sampleStats)
  } else if (contemporaneous == "ggm"){
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta, name = "omega_zeta",
                                                 expcov = contCovEst,
                                                 nNode = nNode,
                                                 nGroup = nGroup,
                                                 labels = vars,
                                                 equal = FALSE, sampletable = sampleStats,
                                                 onlyStartSign = FALSE)
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta, name = "delta_zeta",
                                                 expcov = contCovEst,
                                                 nNode = nNode,
                                                 nGroup = nGroup,
                                                 labels = vars,
                                                 equal = FALSE, sampletable = sampleStats,
                                                 onlyStartSign = FALSE,
                                                 omegaStart = modMatrices$omega_zeta$start)
  } else if (contemporaneous == "prec"){
    modMatrices$kappa_zeta <- matrixsetup_kappa(kappa_zeta, name = "kappa_zeta",
                                                 expcov = contCovEst,
                                                 nNode = nNode,
                                                 nGroup = nGroup,
                                                 labels = vars,
                                                 equal = FALSE, sampletable = sampleStats)
  }

  #### Random effects matrices ####
  # Expected random effects matrix:
  expRanEffects <- as.matrix(spectralshift(sampleStats@covs[[1]] - avgVmat))

  if (randomEffects == "cov"){
    modMatrices$sigma_randomEffects <- matrixsetup_sigma(sigma_randomEffects,
                                                         expcov = list(expRanEffects),
                                                         nNode = nCov,
                                                         nGroup = nGroup,
                                                         labels = covvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "sigma_randomEffects")
  } else if (randomEffects == "chol"){
    modMatrices$lowertri_randomEffects <- matrixsetup_lowertri(lowertri_randomEffects,
                                                               expcov = list(expRanEffects),
                                                               nNode = nCov,
                                                               nGroup = nGroup,
                                                               labels = covvars,
                                                               equal = FALSE,
                                                               sampletable = sampleStats,
                                                               name = "lowertri_randomEffects")
  } else if (randomEffects == "ggm"){
    modMatrices$omega_randomEffects <- matrixsetup_omega(omega_randomEffects,
                                                         expcov = list(expRanEffects),
                                                         nNode = nCov,
                                                         nGroup = nGroup,
                                                         labels = covvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "omega_randomEffects")
    modMatrices$delta_randomEffects <- matrixsetup_delta(delta_randomEffects,
                                                         expcov = list(expRanEffects),
                                                         nNode = nCov,
                                                         nGroup = nGroup,
                                                         labels = covvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "delta_randomEffects",
                                                         omegaStart = modMatrices$omega_randomEffects$start)
  } else if (randomEffects == "prec"){
    modMatrices$kappa_randomEffects <- matrixsetup_kappa(kappa_randomEffects,
                                                         expcov = list(expRanEffects),
                                                         nNode = nCov,
                                                         nGroup = nGroup,
                                                         labels = covvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "kappa_randomEffects")
  } else if (randomEffects == "cor"){
    modMatrices$rho_randomEffects <- matrixsetup_rho(rho_randomEffects,
                                                     expcov = list(expRanEffects),
                                                     nNode = nCov,
                                                     nGroup = nGroup,
                                                     labels = covvars,
                                                     equal = FALSE,
                                                     sampletable = sampleStats,
                                                     name = "rho_randomEffects")
    modMatrices$SD_randomEffects <- matrixsetup_SD(SD_randomEffects,
                                                   expcov = list(expRanEffects),
                                                   nNode = nCov,
                                                   nGroup = nGroup,
                                                   labels = covvars,
                                                   equal = FALSE,
                                                   sampletable = sampleStats,
                                                   name = "SD_randomEffects")
  }

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  # Extra matrices:
  model@extramatrices <- list(
    # VAR1-level matrices (nNode dimension):
    D2 = psychonetrics::duplicationMatrix(nNode),
    L = psychonetrics::eliminationMatrix(nNode),
    Dstar = psychonetrics::duplicationMatrix(nNode, diag = FALSE),
    In = as(diag(nNode), "dMatrix"),
    A = psychonetrics::diagonalizationMatrix(nNode),
    C = as(lavaan::lav_matrix_commutation(nNode, nNode), "dMatrix"),

    # Random effects matrices (nCov dimension):
    D_c = psychonetrics::duplicationMatrix(nCov),
    L_c = psychonetrics::eliminationMatrix(nCov),
    Lstar_c = psychonetrics::eliminationMatrix(nCov, diag = FALSE),
    Dstar_c = psychonetrics::duplicationMatrix(nCov, diag = FALSE),
    In_c = as(diag(nCov), "dMatrix"),
    A_c = psychonetrics::diagonalizationMatrix(nCov),
    C_c = as(lavaan::lav_matrix_commutation(nCov, nCov), "dMatrix"),

    # V matrices:
    V = avgVmat,
    Vall = Vmats,
    Vmethod = Vmethod,
    Vestimation = Vestimation
  )

  # Form the model matrices:
  model@modelmatrices <- formModelMatrices(model)

  ### Baseline model ###
  if (baseline_saturated){
    # Dummy sample stats (not corinput for baseline/saturated):
    sampleStats2 <- sampleStats
    sampleStats2@corinput <- FALSE

    # Form baseline model:
    model@baseline_saturated$baseline <- varcov(metaData,
                                                mu = rep(0, nCov),
                                                type = "chol",
                                                lowertri = "diag",
                                                vars = covvars,
                                                missing = "listwise",
                                                estimator = estimator,
                                                baseline_saturated = FALSE,
                                                sampleStats = sampleStats2)

    model@baseline_saturated$baseline@sample@fullFIML <- FALSE

    ### Saturated model ###
    model@baseline_saturated$saturated <- varcov(metaData,
                                                 type = "chol",
                                                 lowertri = "full",
                                                 vars = covvars,
                                                 missing = "listwise",
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats = sampleStats2)

    model@baseline_saturated$saturated@sample@fullFIML <- FALSE

    # if not FIML, Treat as computed:
    if (estimator != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(
        parVector(model@baseline_saturated$saturated),
        model@baseline_saturated$saturated
      )
    }
  }

  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # Return model:
  return(model)
}

# Wrapper for graphical meta-analytic VAR:
meta_gvar <- function(...){
  meta_var1(..., contemporaneous = "ggm")
}
