# Multi-level lag-1 vector auto-regression (ml_var1).
#
# Two-level (random intercept) VAR(1)/GVAR for intensive longitudinal data,
# estimated by applying the ml_lvm two-level sufficient-statistics ML machinery
# to the lag-embedded ("Toeplitz") data. This is the var1 pseudo-likelihood
# generalized with a between-person level: per subject the lag pairs
# z_ir = (y_{i,t-1}', y_it')' (2p-variate) are treated as conditionally i.i.d.
# given the random intercept b_i, with
#   z_ir | b_i ~ N( mu_z + (b_i',b_i')', Sigma_W ),   b_i ~ N(0, Sigma_B)
#   mu_z    = (mu', mu')'
#   Sigma_W = [[Sigma0, Sigma0 beta'],[beta Sigma0, Sigma0]]   (stationary Toeplitz)
#   Sigma_Bf = 1_{2x2} (x) Sigma_B                              (2p x 2p, rank p)
# and vec(Sigma0) = (I - beta (x) beta)^-1 vec(Sigma_zeta).
#
# The estimand matches panelvar/panelgvar (multilevel VAR with random
# intercepts), but the cost is independent of the number of time points: the
# data are compressed once into the two-level sufficient statistics of the lag
# pairs, and every subsequent evaluation works with 2p x 2p matrices only.
#
# The model layer is implemented in R only (v1); the constructor forces
# usecpp(model, FALSE) so fit/gradient/Fisher route to the R Gauss2L
# estimator (which is model-agnostic and consumes mu, sigma_within,
# sigma_between and the two-level statistics).

# Helper: pooled covariance of the cluster means from complete-data two-level
# sufficient statistics (used only for between-level start values):
.ml_var1_covMeans <- function(ts){
  J <- ts$J
  sizes <- ts$sizes
  p2 <- length(ts$mean_d[[1]])
  grand <- numeric(p2)
  for (s in seq_len(nrow(sizes))) grand <- grand + sizes$m[s] * ts$mean_d[[s]]
  grand <- grand / J
  BB <- matrix(0, p2, p2)
  for (s in seq_len(nrow(sizes))){
    d <- ts$mean_d[[s]] - grand
    BB <- BB + sizes$m[s] * (as.matrix(ts$cov_d[[s]]) + outer(d, d))
  }
  BB / J
}

ml_var1 <- function(
  data,
  idvar,                                  # REQUIRED: subject/cluster variable
  vars,
  dayvar,
  beepvar,
  groups,                                 # deprecated alias, mapped to groupvar
  groupvar,                               # grouping variable
  within_latent  = c("cov","chol","prec","ggm","cor"),
  between_latent = c("cov","chol","prec","ggm","cor"),
  beta = "full",
  omega_zeta_within = "full", delta_zeta_within = "diag",
  kappa_zeta_within = "full", sigma_zeta_within = "full", lowertri_zeta_within = "full",
  rho_zeta_within = "full", SD_zeta_within = "full",
  omega_zeta_between = "full", delta_zeta_between = "diag",
  kappa_zeta_between = "full", sigma_zeta_between = "full", lowertri_zeta_between = "full",
  rho_zeta_between = "full", SD_zeta_between = "full",
  mu,
  equal = "none",
  baseline_saturated = TRUE,
  estimator = "ML",
  optimizer,
  storedata = FALSE,
  standardize = c("none","z","quantile"),
  sampleStats,
  verbose = FALSE,
  toeplitz = TRUE                         # INTERNAL: FALSE = free-within saturated mode
){
  # Experimental flag (fires once per session, unconditional of verbose):
  experimentalWarning("ml_var1 multi-level VAR pseudo-ML estimator")

  # Match args:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  standardize <- match.arg(standardize)

  # CRAN check workarounds:
  . <- NULL

  # Estimator: only "ML" is supported in v1:
  if (!identical(estimator, "ML")){
    stop("ml_var1 only supports estimator = 'ML' (the two-level pseudo-ML estimator). For full-information ML use panelvar()/dlvm1().")
  }

  # groups/groupvar: accept either; standardize to a single 'groups' column name:
  if (missing(groups)) groups <- NULL
  if (!missing(groupvar) && !is.null(groupvar)) groups <- groupvar

  # -----------------------------------------------------------------------
  # Data augmentation + sample statistics (skipped when sampleStats supplied,
  # i.e. for the recursive baseline/saturated calls):
  # -----------------------------------------------------------------------
  if (missing(sampleStats)){

    if (missing(idvar) || is.null(idvar)){
      stop("'idvar' may not be missing for ml_var1 (use var1() for single-level time-series models).")
    }

    if (is.matrix(data)) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("'data' must be a data frame")
    if (is.null(names(data))) names(data) <- paste0("V", seq_len(ncol(data)))

    if (missing(dayvar)) dayvar <- NULL
    if (missing(beepvar)) beepvar <- NULL

    # Resolve variable names:
    if (missing(vars) || is.null(vars)){
      vars <- setdiff(names(data), c(idvar, dayvar, beepvar, groups))
    }

    # Grand (not within-person) standardization of the raw columns before
    # augmentation, so the between variance survives:
    if (standardize == "z"){
      for (v in vars) data[[v]] <- as.numeric(scale(data[[v]], TRUE, TRUE))
    } else if (standardize == "quantile"){
      for (v in vars) data[[v]] <- quantiletransform(data[[v]])
    }

    # Lag-embed (Toeplitz) the data, keeping the subject id per augmented row:
    augData <- tsData(data, vars = vars, beepvar = beepvar, dayvar = dayvar,
                      idvar = idvar, groupvar = groups, includeID = TRUE)

    # Resolve the id column name (tsData uses "ID" when idvar was NULL):
    idcol <- if (!is.null(idvar)) idvar else "ID"

    # The 2p augmented variable columns (lagged first, then current):
    vars2p <- setdiff(colnames(augData), c(groups, idcol))

    # Drop rows with any NA among the 2p variable columns (complete-pairs, v1):
    naRows <- rowSums(is.na(augData[, vars2p, drop = FALSE])) > 0
    if (any(naRows)){
      warning(sum(naRows), " incomplete lag pair(s) with missing values were dropped (complete-pairs analysis).")
      augData <- augData[!naRows, , drop = FALSE]
    }

    # Drop subjects left with zero rows (should not happen after the NA drop,
    # but guard anyway):
    emptySubjects <- setdiff(unique(data[[idvar]]), unique(augData[[idcol]]))
    if (length(emptySubjects) > 0){
      warning(length(emptySubjects), " subject(s) had no complete lag pairs and were dropped.")
    }

    # Sample statistics on the 2p augmented variables (used for bookkeeping and
    # some start values only; estimation uses the two-level statistics below):
    sampleStats <- samplestats(data = augData,
                               vars = vars2p,
                               groups = groups,
                               missing = "listwise",
                               fimldata = FALSE,
                               storedata = storedata,
                               verbose = verbose)

    # Two-level sufficient statistics per group (cluster = subject):
    twolevelStats <- list()
    for (g in seq_len(nrow(sampleStats@groups))){
      if (is.null(groups)){
        subData <- augData
      } else {
        subData <- augData[augData[[groups]] == sampleStats@groups$label[g], , drop = FALSE]
      }
      ts <- twolevel_sufficient_statistics(
        Y = as.matrix(subData[, vars2p, drop = FALSE]),
        cluster = subData[[idcol]]
      )
      if (isTRUE(ts$missing)){
        stop("Internal error: ml_var1 two-level statistics unexpectedly contain missing data after the complete-pairs drop. Please report this bug.")
      }
      twolevelStats[[g]] <- ts
    }
    sampleStats@twolevel <- twolevelStats

    # Override nobs with cluster counts (subjects = independent observations):
    for (g in seq_len(nrow(sampleStats@groups))){
      sampleStats@groups$nobs[g] <- twolevelStats[[g]]$J
    }
  }

  # -----------------------------------------------------------------------
  # Dimensions:
  # -----------------------------------------------------------------------
  nVar2p <- nrow(sampleStats@variables)
  if (nVar2p %% 2 != 0){
    stop("Number of augmented variables is not even; cannot form the Toeplitz structure.")
  }
  p  <- nVar2p / 2                         # number of actual nodes
  kp <- p * (p + 1) / 2                    # vech dim of p x p
  k  <- p * (2 * p + 1)                    # vech dim of 2p x 2p
  labels2p <- sampleStats@variables$label
  labelsNode <- labels2p[p + seq_len(p)]   # current-half labels

  # -----------------------------------------------------------------------
  # Model object:
  # -----------------------------------------------------------------------
  model <- generate_psychonetrics(
    model = "ml_var1",
    submodel = if (within_latent == "ggm" || between_latent == "ggm") "ml_gvar1" else "ml_var1",
    types = list(within_latent = within_latent, between_latent = between_latent, toeplitz = toeplitz),
    sample = sampleStats, computed = FALSE,
    equal = equal,
    optimizer = defaultoptimizer(), estimator = "ML", distribution = "TwoLevelGaussian",
    verbose = verbose)

  nGroup <- nrow(model@sample@groups)

  # Number of "observations" = distribution parameters [mu_z; vech SW; vech SB]
  # per group (so the saturated free-within model has exactly 0 df):
  model@sample@nobs <- nGroup * (2 * p + 2 * k)

  # -----------------------------------------------------------------------
  # Start values (from the two-level pooled-within covariance S_PW):
  # -----------------------------------------------------------------------
  S0est <- list(); betaEst <- list(); zetaWest <- list(); SBest <- list()
  SWfull_start <- list(); SBfull_start <- list()
  twolevelStats <- get_twolevel_stats(sampleStats)
  for (g in seq_len(nGroup)){
    ts <- twolevelStats[[g]]
    S_PW <- as.matrix(ts$S_PW)
    nbar <- ts$N / ts$J

    S0 <- spectralshift(0.5 * (S_PW[1:p, 1:p, drop = FALSE] + S_PW[p + (1:p), p + (1:p), drop = FALSE]))
    S1 <- S_PW[p + (1:p), 1:p, drop = FALSE]
    S0inv <- solve_symmetric(S0)
    bEst <- as.matrix(S1 %*% S0inv)

    S0est[[g]]  <- S0
    betaEst[[g]] <- bEst
    zetaWest[[g]] <- spectralshift(S0 - bEst %*% S0 %*% t(bEst))

    # Between covariance start: pooled cov of cluster means, averaged over the
    # two halves, minus the within-mean sampling variance S0/nbar:
    BB <- .ml_var1_covMeans(ts)
    BBhalf <- 0.5 * (BB[1:p, 1:p, drop = FALSE] + BB[p + (1:p), p + (1:p), drop = FALSE])
    SBest[[g]] <- spectralshift(BBhalf - S0 / nbar)

    # Full 2p start covariances for the saturated (free-within) mode:
    SWfull_start[[g]] <- spectralshift(S_PW)
    SBfull_start[[g]] <- spectralshift(BB - S_PW / nbar)
  }

  # -----------------------------------------------------------------------
  # Model matrices (order fixes the theta ordering!):
  # -----------------------------------------------------------------------
  modMatrices <- list()

  if (toeplitz){
    # mu (p-variate; repeated in both halves by the implied function):
    expMeans <- lapply(model@sample@means, function(x) as.numeric(x)[p + seq_len(p)])
    modMatrices$mu <- matrixsetup_mu(mu, nNode = p, nGroup = nGroup, labels = labelsNode,
                                     equal = "mu" %in% equal, expmeans = expMeans,
                                     sampletable = sampleStats, name = "mu")

    # beta (temporal effects, from = row, to = column as in var1):
    modMatrices$beta <- matrixsetup_beta(beta, name = "beta", nNode = p, nGroup = nGroup,
                                         labels = labelsNode, equal = "beta" %in% equal,
                                         sampletable = sampleStats, start = betaEst,
                                         onlyStartSign = FALSE)

    # Within (contemporaneous) covariance block:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_within, lowertri = lowertri_zeta_within,
                                         omega = omega_zeta_within, delta = delta_zeta_within,
                                         kappa = kappa_zeta_within, rho = rho_zeta_within, SD = SD_zeta_within,
                                         type = within_latent, name = "zeta_within",
                                         sampleStats = sampleStats, equal = equal,
                                         nNode = p, expCov = zetaWest, nGroup = nGroup,
                                         labels = labelsNode))

    # Between (random-intercept) covariance block:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_between, lowertri = lowertri_zeta_between,
                                         omega = omega_zeta_between, delta = delta_zeta_between,
                                         kappa = kappa_zeta_between, rho = rho_zeta_between, SD = SD_zeta_between,
                                         type = between_latent, name = "zeta_between",
                                         sampleStats = sampleStats, equal = equal,
                                         nNode = p, expCov = SBest, nGroup = nGroup,
                                         labels = labelsNode))
  } else {
    # Free-within saturated mode: 2p mu, free 2p within and between blocks (chol):
    expMeans <- lapply(model@sample@means, function(x) as.numeric(x))
    modMatrices$mu <- matrixsetup_mu(mu, nNode = nVar2p, nGroup = nGroup, labels = labels2p,
                                     equal = "mu" %in% equal, expmeans = expMeans,
                                     sampletable = sampleStats, name = "mu")

    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = "full", lowertri = "full",
                                         omega = "full", delta = "diag",
                                         kappa = "full", rho = "full", SD = "full",
                                         type = "chol", name = "zeta_within",
                                         sampleStats = sampleStats, equal = equal,
                                         nNode = nVar2p, expCov = SWfull_start, nGroup = nGroup,
                                         labels = labels2p))

    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = "full", lowertri = "full",
                                         omega = "full", delta = "diag",
                                         kappa = "full", rho = "full", SD = "full",
                                         type = "chol", name = "zeta_between",
                                         sampleStats = sampleStats, equal = equal,
                                         nNode = nVar2p, expCov = SBfull_start, nGroup = nGroup,
                                         labels = labels2p))
  }

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  # -----------------------------------------------------------------------
  # Extra matrices:
  # -----------------------------------------------------------------------
  model@extramatrices <- list(
    # var1 kernel set at dimension p (used by the toeplitz Jacobian):
    In    = as(diag(p), "dMatrix"),
    L     = psychonetrics::eliminationMatrix(p),
    C     = commutationMatrix(p, p),
    D2    = psychonetrics::duplicationMatrix(p),
    Dstar = psychonetrics::duplicationMatrix(p, diag = FALSE),
    A     = psychonetrics::diagonalizationMatrix(p),

    # Estimator-layer duplication matrix at dimension 2p (consumed by
    # expected_hessian_Gauss2L_group):
    D_y   = psychonetrics::duplicationMatrix(nVar2p),

    # 2p-dimension kernels (used by the saturated free-within Jacobian):
    I2p   = as(diag(nVar2p), "dMatrix"),
    L2p   = psychonetrics::eliminationMatrix(nVar2p),
    C2p   = commutationMatrix(nVar2p, nVar2p)
  )

  # Structure maps (dummy-index construction, one 1 per row):
  # P_within: vech(Sigma_W) = P_within %*% [ vech(Sigma0) ; vec(Sigma1) ]
  d0 <- matrix(0, p, p); d0[lower.tri(d0, diag = TRUE)] <- seq_len(kp)
  d0[upper.tri(d0)] <- t(d0)[upper.tri(d0)]
  d1 <- matrix(kp + seq_len(p^2), p, p)
  SWdummy <- rbind(cbind(d0, matrix(0, p, p)), cbind(d1, d0))
  inds <- SWdummy[lower.tri(SWdummy, diag = TRUE)]
  model@extramatrices$P_within <- as(
    sparseMatrix(i = seq_along(inds), j = inds, dims = c(k, kp + p^2)), "dMatrix")

  # P_between: vech(1_{2x2} (x) Sigma_B) = P_between %*% vech(Sigma_B)
  SBdummy <- matrix(1, 2, 2) %x% d0
  indsB <- SBdummy[lower.tri(SBdummy, diag = TRUE)]
  model@extramatrices$P_between <- as(
    sparseMatrix(i = seq_along(indsB), j = indsB, dims = c(k, kp)), "dMatrix")

  # Form the model matrices:
  model@modelmatrices <- formModelMatrices(model)

  # -----------------------------------------------------------------------
  # Baseline / saturated models:
  # -----------------------------------------------------------------------
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){

    # Baseline (independence): beta = 0, diagonal within and between blocks:
    model@baseline_saturated$baseline <- ml_var1(
      idvar = idvar,
      beta = "zero",
      within_latent = "chol", lowertri_zeta_within = "diag",
      between_latent = "chol", lowertri_zeta_between = "diag",
      equal = equal,
      estimator = "ML",
      baseline_saturated = FALSE,
      sampleStats = sampleStats,
      toeplitz = TRUE)

    # Saturated: free-within mode (no equal; unconstrained per group). The
    # free 2p within and between blocks use a Cholesky parameterization, so the
    # types must be "chol" (the derivatives/implied dispatch on them):
    model@baseline_saturated$saturated <- ml_var1(
      idvar = idvar,
      within_latent = "chol", between_latent = "chol",
      estimator = "ML",
      baseline_saturated = FALSE,
      sampleStats = sampleStats,
      toeplitz = FALSE)

    # For estimator "ML" the saturated model has no closed form: leave it
    # computed = FALSE so runmodel() optimizes it (mirrors ml_lvm).
  }

  # Identify (no-op for ml_var1):
  model <- identify(model)

  # Optimizer:
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # v1 implements the model layer in R only: force the R Gauss2L path:
  model <- usecpp(model, FALSE)

  return(model)
}

# GGM-defaults wrapper (mirrors gvar/panelgvar):
ml_gvar1 <- function(..., within_latent = "ggm", between_latent = "ggm"){
  ml_var1(..., within_latent = within_latent, between_latent = between_latent)
}
