# Multi-level variance-covariance family (ml_varcov) and its GGM / correlation
# wrappers (ml_ggm / ml_corr, in their own files).
#
# ml_varcov is the multi-level analogue of the single-level varcov family: it
# models the within-cluster and between-cluster variance-covariance matrices of
# clustered (two-level) data directly, with no latent layer. Each level can be
# parameterized independently as a covariance matrix ("cov"), a Cholesky
# decomposition ("chol"), a precision matrix ("prec"), a Gaussian graphical
# model ("ggm"; the multi-level GGM) or a correlation matrix ("cor"), exactly as
# the single-level varcov 'type' argument does.
#
# The model is estimated with the two-level sufficient-statistics Gaussian ML
# estimator (distribution "TwoLevelGaussian"), the same estimator ml_lvm /
# ml_var1 use: the cost is independent of the number of units per cluster. The
# model layer (implied + derivatives) is implemented natively for ml_varcov and
# carries only mu, sigma_within and sigma_between -- there is no lambda, beta,
# residual or latent-mean parameter. The single-level varcov covariance-
# structure derivatives (d_sigma_omega, d_sigma_delta, ...) are reused per level.

# Helper: pooled covariance of the cluster means from complete-data two-level
# sufficient statistics (used only for the between-level start values):
.ml_varcov_covMeans <- function(ts){
  J <- ts$J
  sizes <- ts$sizes
  p <- length(ts$mean_d[[1]])
  grand <- numeric(p)
  for (s in seq_len(nrow(sizes))) grand <- grand + sizes$m[s] * ts$mean_d[[s]]
  grand <- grand / J
  BB <- matrix(0, p, p)
  for (s in seq_len(nrow(sizes))){
    d <- ts$mean_d[[s]] - grand
    BB <- BB + sizes$m[s] * (as.matrix(ts$cov_d[[s]]) + outer(d, d))
  }
  BB / J
}

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
  estimator = "ML",
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

  # CRAN check workaround:
  . <- NULL

  # Estimator: only the two-level ML estimator is supported:
  if (!identical(estimator, "ML")){
    stop("ml_varcov only supports estimator = 'ML' (the two-level sufficient-statistics ML estimator).")
  }

  # groups/groupvar: accept either; standardize to a single 'groups' column name:
  if (missing(groups)) groups <- NULL
  if (!missing(groupvar) && !is.null(groupvar)) groups <- groupvar

  # -----------------------------------------------------------------------
  # Data + sample statistics (skipped when sampleStats supplied, i.e. for the
  # recursive baseline/saturated calls):
  # -----------------------------------------------------------------------
  if (missing(sampleStats)){

    if (missing(clusters) || is.null(clusters)){
      stop("'clusters' may not be missing for ml_varcov (use the varcov() family for single-level models).")
    }

    if (is.matrix(data)) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("'data' must be a data frame")
    if (is.null(names(data))) stop("Dataset contains no column names.")

    # Resolve variable names:
    if (missing(vars) || is.null(vars)){
      vars <- setdiff(names(data), c(clusters, groups))
    }

    if (!clusters %in% names(data)){
      stop("'clusters' argument does not correspond to a column name of 'data'")
    }

    # Remove rows with NA cluster:
    if (any(is.na(data[[clusters]]))){
      warning("Rows with NA cluster removed.")
      data <- data[!is.na(data[[clusters]]), , drop = FALSE]
    }

    # Grand (not within-cluster) standardization of the raw columns before
    # forming statistics, so the between variance survives:
    if (standardize == "z"){
      for (v in vars) data[[v]] <- as.numeric(scale(data[[v]], TRUE, TRUE))
    } else if (standardize == "quantile"){
      for (v in vars) data[[v]] <- quantiletransform(data[[v]])
    }

    # Grouping column (single-group default):
    if (is.null(groups)){
      groups <- "GROUPID"
      data[[groups]] <- "fullsample"
    }

    # Sample statistics on the raw long data (used for group/variable
    # bookkeeping and start values only; estimation uses the two-level
    # statistics below). Muffle the harmless NA-covariance warning that very
    # small groups can trigger:
    sampleStats <- withCallingHandlers(
      samplestats(data = data,
                  vars = vars,
                  groups = groups,
                  missing = "listwise",
                  fimldata = FALSE,
                  storedata = storedata,
                  verbose = verbose),
      warning = function(w){
        if (grepl("NA sample covariances", conditionMessage(w))){
          invokeRestart("muffleWarning")
        }
      })

    # Two-level sufficient statistics per group (cluster = clusters):
    twolevelStats <- list()
    for (g in seq_len(nrow(sampleStats@groups))){
      subData <- data[data[[groups]] == sampleStats@groups$label[g], , drop = FALSE]
      ts <- twolevel_sufficient_statistics(
        Y = as.matrix(subData[, vars, drop = FALSE]),
        cluster = subData[[clusters]]
      )
      if (isTRUE(ts$missing)){
        stop("ml_varcov does not yet support missing data; please supply complete cases.")
      }
      twolevelStats[[g]] <- ts
    }
    sampleStats@twolevel <- twolevelStats

    # The number of independent observations equals the number of clusters:
    for (g in seq_len(nrow(sampleStats@groups))){
      sampleStats@groups$nobs[g] <- twolevelStats[[g]]$J
    }
  }

  # -----------------------------------------------------------------------
  # Dimensions:
  # -----------------------------------------------------------------------
  nVar <- nrow(sampleStats@variables)
  kp <- nVar * (nVar + 1) / 2
  labels <- sampleStats@variables$label

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
    optimizer = defaultoptimizer(), estimator = "ML", distribution = "TwoLevelGaussian",
    verbose = verbose)

  nGroup <- nrow(model@sample@groups)

  # Number of "observations" = distribution parameters [mu; vech SW; vech SB]
  # per group (so the saturated model has exactly 0 df):
  model@sample@nobs <- nGroup * (nVar + 2 * kp)

  # -----------------------------------------------------------------------
  # Start values (within = pooled-within covariance; between = pooled cov of
  # cluster means minus the within sampling variance):
  # -----------------------------------------------------------------------
  withinStart <- list(); betweenStart <- list()
  twolevelStats <- get_twolevel_stats(sampleStats)
  for (g in seq_len(nGroup)){
    ts <- twolevelStats[[g]]
    S_PW <- spectralshift(as.matrix(ts$S_PW))
    nbar <- ts$N / ts$J
    BB <- .ml_varcov_covMeans(ts)
    withinStart[[g]]  <- S_PW
    betweenStart[[g]] <- spectralshift(BB - S_PW / nbar)
  }

  # -----------------------------------------------------------------------
  # Model matrices (order fixes the theta ordering!):
  # -----------------------------------------------------------------------
  modMatrices <- list()

  # mu (p-variate):
  expMeans <- lapply(model@sample@means, function(x) as.numeric(x))
  modMatrices$mu <- matrixsetup_mu(mu, nNode = nVar, nGroup = nGroup, labels = labels,
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
                                       labels = labels))

  # Between-cluster covariance block:
  modMatrices <- c(modMatrices,
                   matrixsetup_flexcov(sigma = sigma_between, lowertri = lowertri_between,
                                       omega = omega_between, delta = delta_between,
                                       kappa = kappa_between, rho = rho_between, SD = SD_between,
                                       type = between, name = "between",
                                       sampleStats = sampleStats, equal = equal,
                                       nNode = nVar, expCov = betweenStart, nGroup = nGroup,
                                       labels = labels))

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

  # Form the model matrices:
  model@modelmatrices <- formModelMatrices(model)

  # -----------------------------------------------------------------------
  # Baseline / saturated models:
  # -----------------------------------------------------------------------
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){

    # Baseline (independence): diagonal within and between covariance blocks:
    model@baseline_saturated$baseline <- ml_varcov(
      within = "chol", lowertri_within = "diag",
      between = "chol", lowertri_between = "diag",
      equal = equal,
      estimator = "ML",
      baseline_saturated = FALSE,
      sampleStats = sampleStats)

    # Saturated: free within and between covariance blocks (Cholesky
    # parameterization; the derivatives/implied dispatch on the types):
    model@baseline_saturated$saturated <- ml_varcov(
      within = "chol", between = "chol",
      estimator = "ML",
      baseline_saturated = FALSE,
      sampleStats = sampleStats)

    # For estimator "ML" the saturated model has no closed form: leave it
    # computed = FALSE so runmodel() optimizes it (mirrors ml_lvm / ml_var1).
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
  # against the R path). usecpp() selects the C++ path where available:
  model <- usecpp(model)

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
