# Multi-level variance-covariance family (ml_varcov) and its GGM / correlation
# wrappers (ml_ggm / ml_corr).
#
# ml_varcov is the multi-level analogue of the single-level varcov family: it
# models the within-cluster and between-cluster variance-covariance matrices of
# clustered (two-level) data directly, with no latent layer. Each level can be
# parameterized independently as a covariance matrix ("cov"), a Cholesky
# decomposition ("chol"), a precision matrix ("prec"), a Gaussian graphical
# model ("ggm"; the multi-level GGM) or a correlation matrix ("cor"), exactly as
# the single-level varcov 'type' argument does.
#
# Two estimators are supported, mirroring ml_lvm:
# - "ML": the two-level sufficient-statistics Gaussian ML estimator
#   (distribution "TwoLevelGaussian"); the cost is independent of the number of
#   units per cluster.
# - "FIML": wide-format full-information ML (distribution "Gaussian" with
#   per-pattern FIML data), the same likelihood evaluated on the one-row-per-
#   cluster wide format. Supports missing data through the FIML patterns.
# The model layer (implied + derivatives) is implemented natively for ml_varcov
# and carries only mu, sigma_within and sigma_between -- there is no lambda,
# beta, residual or latent-mean parameter. The single-level varcov covariance-
# structure derivatives (d_sigma_omega, d_sigma_delta, ...) are reused per level.

ml_varcov <- function(
  data,
  clusters,                                 # REQUIRED: cluster variable
  vars,
  within  = c("cov","chol","prec","ggm","cor"),
  between = c("cov","chol","prec","ggm","cor"),
  groups,                                   # deprecated alias, mapped to groupvar
  groupvar,                                 # grouping variable
  # Within-cluster structure specifications:
  sigma_within = "full", lowertri_within = "full", omega_within = "full",
  delta_within = "diag", kappa_within = "full", rho_within = "full", SD_within = "full",
  # Between-cluster structure specifications:
  sigma_between = "full", lowertri_between = "full", omega_between = "full",
  delta_between = "diag", kappa_between = "full", rho_between = "full", SD_between = "full",
  mu,
  equal = "none",
  baseline_saturated = TRUE,
  estimator = c("default","ML","FIML"),
  optimizer,
  storedata = FALSE,
  standardize = c("none","z","quantile"),
  sampleStats,
  verbose = FALSE
){
  # Experimental flag (fires once per session, unconditional of verbose):
  experimentalWarning("ml_varcov multi-level variance-covariance / GGM model")

  # Match args:
  within  <- match.arg(within)
  between <- match.arg(between)
  standardize <- match.arg(standardize)
  estimator <- match.arg(estimator)

  # CRAN check workaround:
  . <- NULL

  # groups/groupvar: accept either; standardize to a single 'groups' column name:
  if (missing(groups)) groups <- NULL
  if (!missing(groupvar) && !is.null(groupvar)) groups <- groupvar

  # -----------------------------------------------------------------------
  # Data processing (always; the recursive baseline/saturated calls pass the
  # already-processed data plus 'sampleStats', in which case only the sample-
  # statistics computation below is skipped):
  # -----------------------------------------------------------------------
  if (missing(clusters) || is.null(clusters)){
    stop("'clusters' may not be missing for ml_varcov (use the varcov() family for single-level models).")
  }

  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (is.null(names(data))) stop("Dataset contains no column names.")

  # Resolve variable names:
  if (missing(vars) || is.null(vars)){
    vars <- setdiff(names(data), c(clusters, groups, "CLUSTERID"))
  }
  nVar <- length(vars)
  kp <- nVar * (nVar + 1) / 2

  if (!clusters %in% names(data)){
    stop("'clusters' argument does not correspond to a column name of 'data'")
  }

  # Remove rows with NA cluster:
  if (any(is.na(data[[clusters]]))){
    warning("Rows with NA cluster removed.")
    data <- data[!is.na(data[[clusters]]), , drop = FALSE]
  }

  # Grouping column (single-group default; resolved before the cluster IDs so
  # that cluster labels reused across groups denote different clusters):
  if (is.null(groups)){
    groups <- "GROUPID"
    data[[groups]] <- "fullsample"
  }
  groupLabels <- unique(data[[groups]])

  # Add a column with ID in cluster (used by the wide format). ave() assigns
  # in the original row order, so the data need not be sorted by cluster:
  data[['CLUSTERID']] <- ave(seq_len(nrow(data)), data[[groups]], data[[clusters]], FUN = seq_along)
  maxInCluster <- max(data[['CLUSTERID']])
  anyMissing <- anyNA(data[, vars])

  # Resolve the estimator (mirrors ml_lvm):
  if (estimator == "default"){
    if (!anyMissing && maxInCluster > 5){
      estimator <- "ML"
    } else {
      estimator <- "FIML"
    }
    if (verbose){
      message("Using estimator = '", estimator, "' (set the 'estimator' argument to overwrite).")
    }
  }
  if (!estimator %in% c("ML","FIML")){
    stop("ml_varcov only supports estimator = 'ML' (two-level sufficient-statistics ML) or estimator = 'FIML' (wide-format full-information ML).")
  }
  if (estimator == "ML" && anyMissing && verbose){
    experimentalWarning("ml_varcov two-level ML estimator with within-cluster missing data")
  }

  # Grand (not within-cluster) standardization of the raw columns before
  # forming statistics, so the between variance survives. The recursive
  # baseline/saturated calls receive the already-standardized data and the
  # default standardize = "none", so no double transformation occurs:
  if (standardize == "z"){
    for (v in vars) data[[v]] <- as.numeric(scale(data[[v]], TRUE, TRUE))
  } else if (standardize == "quantile"){
    for (v in vars) data[[v]] <- quantiletransform(data[[v]])
  }

  # -----------------------------------------------------------------------
  # Start values per group (missing-data proof): the within start is the
  # covariance of the cluster-centered data, the between start the covariance
  # of the cluster means minus the within-mean sampling variance:
  # -----------------------------------------------------------------------
  withinStart <- list(); betweenStart <- list()
  for (g in seq_along(groupLabels)){
    subData <- data[data[[groups]] == groupLabels[g], , drop = FALSE]
    Y <- as.matrix(subData[, vars, drop = FALSE])
    cl <- as.character(subData[[clusters]])
    cm <- do.call(cbind, lapply(seq_len(nVar), function(i) tapply(Y[, i], cl, mean, na.rm = TRUE)))
    Yc <- Y - cm[match(cl, rownames(cm)), , drop = FALSE]
    SW_s <- spectralshift(cov(Yc, use = "pairwise.complete.obs"))
    nbar <- nrow(Y) / nrow(cm)
    SB_s <- spectralshift(cov(cm, use = "pairwise.complete.obs") - SW_s / nbar)
    withinStart[[g]] <- SW_s
    betweenStart[[g]] <- SB_s
  }

  # -----------------------------------------------------------------------
  # Wide format (FIML only): one row per cluster, columns var_position:
  # -----------------------------------------------------------------------
  if (estimator == "FIML"){
    datalong <- tidyr::gather(data[, c(vars, clusters, groups, "CLUSTERID")], "variable", "value", vars)
    datawide <- tidyr::pivot_wider(datalong, id_cols = c(clusters, groups),
                                   values_from = "value", names_from = c("variable", "CLUSTERID"))

    # Design matrix (exact matching: variable names may contain regex
    # metacharacters, mirroring ml_lvm):
    design <- matrix(NA_character_, nVar, maxInCluster)
    for (i in seq_len(nVar)){
      for (j in seq_len(maxInCluster)){
        varName <- paste0(vars[i], "_", j)
        if (sum(names(datawide) == varName) == 1){
          design[i, j] <- varName
        }
      }
    }
    allVars <- na.omit(as.vector(design))
    nAllVar <- length(allVars)
    designPattern <- as(1 * (!is.na(design)), "matrix")
    casesPerVar <- as.vector(designPattern * row(designPattern))
    casesPerVar <- casesPerVar[casesPerVar != 0]
  }

  # -----------------------------------------------------------------------
  # Sample statistics (skipped for the recursive baseline/saturated calls):
  # -----------------------------------------------------------------------
  if (missing(sampleStats)){
    if (estimator == "FIML"){
      # Wide-format sample statistics with per-pattern FIML data:
      sampleStats <- samplestats(data = datawide,
                                 vars = allVars,
                                 groups = groups,
                                 fimldata = TRUE,
                                 storedata = storedata,
                                 verbose = verbose)
    } else {
      # Long-format sample statistics (bookkeeping and start values only;
      # estimation uses the two-level sufficient statistics below). Muffle the
      # harmless NA-covariance warning that small groups can trigger:
      sampleStats <- withCallingHandlers(
        samplestats(data = data,
                    vars = vars,
                    groups = groups,
                    missing = "pairwise",
                    fimldata = FALSE,
                    storedata = storedata,
                    verbose = verbose),
        warning = function(w){
          if (grepl("NA sample covariances", conditionMessage(w))){
            invokeRestart("muffleWarning")
          }
        })

      # Two-level sufficient statistics per group:
      twolevelStats <- list()
      for (g in seq_len(nrow(sampleStats@groups))){
        subData <- data[data[[groups]] == sampleStats@groups$label[g], , drop = FALSE]
        twolevelStats[[g]] <- twolevel_sufficient_statistics(
          Y = as.matrix(subData[, vars, drop = FALSE]),
          cluster = subData[[clusters]]
        )
      }
      sampleStats@twolevel <- twolevelStats

      # The number of independent observations equals the number of clusters:
      for (g in seq_len(nrow(sampleStats@groups))){
        sampleStats@groups$nobs[g] <- twolevelStats[[g]]$J
      }
    }
  }

  # -----------------------------------------------------------------------
  # Model object:
  # -----------------------------------------------------------------------
  submodel <- if (within == "ggm" && between == "ggm"){
    "ml_ggm"
  } else if (within == "cor" && between == "cor"){
    "ml_corr"
  } else {
    "ml_varcov"
  }

  model <- generate_psychonetrics(
    model = "ml_varcov",
    submodel = submodel,
    types = list(within = within, between = between),
    sample = sampleStats, computed = FALSE,
    equal = equal,
    optimizer = defaultoptimizer(), estimator = estimator,
    distribution = ifelse(estimator == "ML", "TwoLevelGaussian", "Gaussian"),
    verbose = verbose)

  nGroup <- nrow(model@sample@groups)

  # Number of "observations":
  if (estimator == "ML"){
    # Distribution parameters [mu; vech SW; vech SB] per group (so the
    # saturated model has exactly 0 df):
    model@sample@nobs <- nGroup * (nVar + 2 * kp)
  } else {
    # Wide-format means and covariances per group (mirrors ml_lvm; the wide
    # saturated correction cancels in the df computation):
    model@sample@nobs <- nAllVar * (nAllVar + 1) / 2 * nGroup + nAllVar * nGroup
  }

  # -----------------------------------------------------------------------
  # Model matrices (order fixes the theta ordering!):
  # -----------------------------------------------------------------------
  modMatrices <- list()

  # Expected means (per group, p-variate):
  if (estimator == "FIML"){
    expMeans <- lapply(model@sample@means, function(x) as.numeric(tapply(as.numeric(x), casesPerVar, mean, na.rm = TRUE)))
  } else {
    expMeans <- lapply(model@sample@means, function(x) as.numeric(x))
  }

  # mu (p-variate):
  modMatrices$mu <- matrixsetup_mu(mu, nNode = nVar, nGroup = nGroup, labels = vars,
                                   equal = "mu" %in% equal, expmeans = expMeans,
                                   sampletable = sampleStats, name = "mu")

  # Within-cluster covariance block:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma = sigma_within, lowertri = lowertri_within,
                                       omega = omega_within, delta = delta_within,
                                       kappa = kappa_within, rho = rho_within, SD = SD_within,
                                       type = within, name = "within",
                                       sampleStats = sampleStats, equal = equal,
                                       nNode = nVar, expCov = withinStart, nGroup = nGroup,
                                       labels = vars))

  # Between-cluster covariance block:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma = sigma_between, lowertri = lowertri_between,
                                       omega = omega_between, delta = delta_between,
                                       kappa = kappa_between, rho = rho_between, SD = SD_between,
                                       type = between, name = "between",
                                       sampleStats = sampleStats, equal = equal,
                                       nNode = nVar, expCov = betweenStart, nGroup = nGroup,
                                       labels = vars))

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  # -----------------------------------------------------------------------
  # Extra matrices (p-dimensional derivative kernels + estimator duplication):
  # -----------------------------------------------------------------------
  model@extramatrices <- list(
    In    = as(diag(nVar), "dMatrix"),
    L     = psychonetrics::eliminationMatrix(nVar),
    C     = commutationMatrix(nVar, nVar),
    D     = psychonetrics::duplicationMatrix(nVar),
    Dstar = psychonetrics::duplicationMatrix(nVar, diag = FALSE),
    A     = psychonetrics::diagonalizationMatrix(nVar),
    # Consumed by expected_hessian_Gauss2L_group:
    D_y   = psychonetrics::duplicationMatrix(nVar)
  )

  if (estimator == "FIML"){
    # Design pattern and the permutation P mapping the compact model rows
    # [mu (p); vech within-position block (kp); vec between-positions block
    # (p^2)] to the wide-format distribution parameters [wide mu; vech wide
    # Sigma] (mirrors the ml_lvm construction):
    model@extramatrices$designPattern <- designPattern

    muDummy <- matrix(rep(seq_len(nVar), maxInCluster))
    sigDummy <- matrix(0, nVar, nVar)
    sigDummy[lower.tri(sigDummy, diag = TRUE)] <- max(muDummy) + seq_len(kp)
    sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]

    U <- list(sigDummy)
    if (maxInCluster > 1){
      U <- c(U, lapply(seq_len(maxInCluster - 1), function(x) matrix(max(sigDummy) + seq_len(nVar^2), nVar, nVar)))
    }
    allSigmas <- blockToeplitz(U)
    totElements <- max(allSigmas)

    subMu <- muDummy[as.vector(designPattern == 1), , drop = FALSE]
    subSigmas <- allSigmas[as.vector(designPattern == 1), as.vector(designPattern == 1)]

    distVec <- c(as.vector(subMu), subSigmas[lower.tri(subSigmas, diag = TRUE)])
    nTotal <- length(distVec)
    distVecrawts <- seq_along(distVec)[distVec != 0]
    distVec <- distVec[distVec != 0]

    model@extramatrices$P <- as(sparseMatrix(
      i = distVecrawts, j = distVec, dims = c(nTotal, totElements)
    ), "dMatrix")
  }

  # Form the model matrices:
  model@modelmatrices <- formModelMatrices(model)

  # -----------------------------------------------------------------------
  # Baseline / saturated models (recursive calls; the processed data are
  # passed along with the sample statistics):
  # -----------------------------------------------------------------------
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){

    # Baseline (independence): diagonal within and between covariance blocks:
    model@baseline_saturated$baseline <- ml_varcov(
      data = data, clusters = clusters, vars = vars, groupvar = groups,
      within = "chol", lowertri_within = "diag",
      between = "chol", lowertri_between = "diag",
      equal = equal,
      estimator = estimator,
      baseline_saturated = FALSE,
      sampleStats = sampleStats)

    # Saturated: free within and between covariance blocks (Cholesky
    # parameterization; the derivatives/implied dispatch on the types). No
    # 'equal': the saturated reference is unconstrained per group:
    model@baseline_saturated$saturated <- ml_varcov(
      data = data, clusters = clusters, vars = vars, groupvar = groups,
      within = "chol", between = "chol",
      estimator = estimator,
      baseline_saturated = FALSE,
      sampleStats = sampleStats)

    # For estimators "ML" and "FIML" the saturated model has no closed form:
    # leave it computed = FALSE so runmodel() optimizes it (mirrors ml_lvm).
  }

  # Identify (no-op for ml_varcov -- no latent scales to fix):
  model <- identify(model)

  # Optimizer:
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # The model layer is implemented in both R and C++ (the C++ twin is verified
  # against the R path); usecpp() selects the C++ path where available. The
  # two-level ML estimator with within-cluster missing data uses the R path
  # (the per-pattern missing-data likelihood has no C++ twin) and numeric
  # Fisher information, mirroring ml_lvm:
  if (estimator == "ML" && anyMissing){
    model <- usecpp(model, FALSE)
  } else {
    model <- usecpp(model)
  }

  return(model)
}

# GGM-defaults wrapper: multi-level Gaussian graphical model (both levels GGM).
ml_ggm <- function(...){
  ml_varcov(..., within = "ggm", between = "ggm")
}

# Correlation-defaults wrapper: multi-level correlation model (both levels cor).
ml_corr <- function(...){
  ml_varcov(..., within = "cor", between = "cor")
}
