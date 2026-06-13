# Add fit measures to psychonetrics object!


# Helper to safely check if WLS.Gamma slot exists and is non-empty:
has_WLS_Gamma <- function(x) {
  .hasSlot(x@sample, "WLS.Gamma") && length(x@sample@WLS.Gamma) > 0
}

# Build the full model Jacobian Delta = d_phi_theta_<model>(prep) %*% Mmatrix.
# This is the derivative of the (stacked, per-group) sample-statistic vector
# (means then vech(Sigma), column-major lower triangle) with respect to the
# free parameters. Shared by the least-squares sandwich (wls_sandwich_components)
# and the robust ML machinery (ml_robust_components).
build_Delta_full <- function(x){
  if (x@cpp){
    prep <- prepareModel_cpp(parVector(x), x)
  } else {
    prep <- prepareModel(parVector(x), x)
  }
  if (x@cpp){
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov_cpp,
      "lvm" = d_phi_theta_lvm_cpp,
      "var1" = d_phi_theta_var1_cpp,
      "dlvm1" = d_phi_theta_dlvm1_cpp,
      "tsdlvm1" = d_phi_theta_tsdlvm1_cpp,
      "meta_varcov" = d_phi_theta_meta_varcov_cpp,
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising_cpp,
      "ml_lvm" = d_phi_theta_ml_lvm_cpp,
      "meta_lvm" = d_phi_theta_meta_lvm_cpp,
      "meta_var1" = d_phi_theta_meta_var1_cpp
    )
  } else {
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov,
      "lvm" = d_phi_theta_lvm,
      "var1" = d_phi_theta_var1,
      "dlvm1" = d_phi_theta_dlvm1,
      "tsdlvm1" = d_phi_theta_tsdlvm1,
      "meta_varcov" = d_phi_theta_meta_varcov,
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising,
      "ml_lvm" = d_phi_theta_ml_lvm,
      "meta_lvm" = d_phi_theta_meta_lvm,
      "meta_var1" = d_phi_theta_meta_var1
    )
  }
  modelPart <- modelJacobian(prep)
  if (x@cpp){
    M <- Mmatrix_cpp(x@parameters)
  } else {
    M <- Mmatrix(x@parameters)
  }
  as.matrix(modelPart %*% M)
}

# Normal-theory (ML) weight matrix V for one group, in the sample-statistic
# ordering (means then vech(Sigma)). V = bdiag(Sigma^-1, 0.5 D'(Sigma^-1 (x)
# Sigma^-1) D). Honours meanstructure / corinput row dropping exactly as
# LS_weightsmat() does (drop the leading nvar mean rows when !meanstructure;
# drop the variance rows when corinput). When kappa = Sigma^-1 is already
# available it is reused to avoid a re-inversion.
ml_Vmat <- function(Sigma, meanstructure = TRUE, corinput = FALSE, kappa = NULL){
  Sigma <- as.matrix(Sigma)
  p <- ncol(Sigma)
  Si <- if (is.null(kappa)) solve(Sigma) else as.matrix(kappa)
  D <- lavaan::lav_matrix_duplication(p)
  V22 <- as.matrix(0.5 * t(D) %*% (Si %x% Si) %*% D)
  nstat <- p + p * (p + 1) / 2
  V <- matrix(0, nstat, nstat)
  V[seq_len(p), seq_len(p)] <- Si
  V[-seq_len(p), -seq_len(p)] <- V22
  # Drop mean rows/cols if no mean structure:
  if (!meanstructure){
    V <- V[-seq_len(p), -seq_len(p), drop = FALSE]
  }
  # Drop variance rows/cols if correlation input:
  if (corinput){
    inds <- meanstructure * p + which(diag(p)[lower.tri(diag(p), diag = TRUE)] == 1)
    V <- V[-inds, -inds, drop = FALSE]
  }
  V
}

# Assemble the per-group robust-ML building blocks for estimator = "ML":
#   Delta_g : group rows of the model Jacobian (means then vech(Sigma))
#   V_g     : normal-theory weight matrix from the model-implied Sigma_g
#   Gamma_g : asymptotic covariance of the sample statistics (x@sample@WLS.Gamma)
#   E       : x@information (unit expected Fisher information; SEs from
#             solve(E)/N reproduce lavaan se = "standard")
#   fg      : group weights n_g / N
# Returns NULL when the estimator is not ML or Gamma is unavailable.
ml_robust_components <- function(x){
  if (x@estimator != "ML") return(NULL)
  if (!has_WLS_Gamma(x)) return(NULL)
  if (length(x@modelmatrices) == 0) return(NULL)

  nGroups <- nrow(x@sample@groups)
  nobs_per_group <- x@sample@groups$nobs
  nTotal <- sum(nobs_per_group)

  Delta_full <- build_Delta_full(x)

  meanstructure <- x@meanstructure
  corinput <- isTRUE(x@sample@corinput)

  # Row count per group from the stored Gamma:
  nrows_per_group <- sapply(x@sample@WLS.Gamma, nrow)
  if (sum(nrows_per_group) != nrow(Delta_full)) return(NULL)

  Delta_list <- vector("list", nGroups)
  V_list <- vector("list", nGroups)
  Gamma_list <- vector("list", nGroups)
  row_offset <- 0
  for (g in seq_len(nGroups)){
    ng <- nrows_per_group[g]
    row_inds <- row_offset + seq_len(ng)
    Delta_list[[g]] <- Delta_full[row_inds, , drop = FALSE]
    row_offset <- row_offset + ng

    Sigma_g <- x@modelmatrices[[g]]$sigma
    kappa_g <- x@modelmatrices[[g]]$kappa
    V_list[[g]] <- ml_Vmat(Sigma_g, meanstructure = meanstructure,
                           corinput = corinput, kappa = kappa_g)
    Gamma_list[[g]] <- as.matrix(x@sample@WLS.Gamma[[g]])
  }

  # Expected unit information E (= x@information). Recompute if absent.
  if (!is.null(x@information) && length(x@information) > 0){
    E <- as.matrix(x@information)
  } else if (x@cpp){
    E <- psychonetrics_FisherInformation_cpp(x)
  } else {
    E <- psychonetrics_FisherInformation(x)
  }

  list(
    Delta = Delta_list,
    V = V_list,
    Gamma = Gamma_list,
    E = E,
    fg = nobs_per_group / nTotal,
    nGroups = nGroups,
    npar = ncol(Delta_full),
    nTotal = nTotal
  )
}

# Casewise scores of the saturated (h1) Gaussian log-likelihood w.r.t. the
# sample statistics (mu, vech(Sigma)), evaluated at the supplied mu/Sigma.
# Rows = cases, columns = (means then vech(Sigma)) in psychonetrics ordering.
# Honours meanstructure / corinput row dropping exactly as ml_Vmat() / LS_Gamma.
#   mean part:  w_i = Sigma^-1 (y_i - mu)
#   cov  part:  m_i,(a,b) = 0.5 (w_ia w_ib - (Sigma^-1)_ab), with the duplication
#               (vech) factor applied by halving the diagonal entries.
ml_casewise_scores_h1 <- function(Y, mu, Sigma, meanstructure = TRUE, corinput = FALSE){
  Y <- as.matrix(Y)
  p <- ncol(Y)
  Si <- solve(Sigma)
  Yc <- sweep(Y, 2, mu)            # y_i - mu
  W <- Yc %*% Si                   # dlogl_i/dmu = Si (y_i - mu)
  ii <- unlist(lapply(seq_len(p), function(j) j:p))   # vech row idx (col-major)
  jj <- rep(seq_len(p), times = p:1)                  # vech col idx
  Z <- W[, ii, drop = FALSE] * W[, jj, drop = FALSE]
  iS <- Si[cbind(ii, jj)]
  Sc_s <- sweep(Z, 2, iS)
  diagh <- which(ii == jj)
  Sc_s[, diagh] <- Sc_s[, diagh] / 2
  SC <- cbind(W, Sc_s)
  if (!meanstructure){
    SC <- SC[, -seq_len(p), drop = FALSE]
  }
  if (corinput){
    inds <- meanstructure * p + which(diag(p)[lower.tri(diag(p), diag = TRUE)] == 1)
    SC <- SC[, -inds, drop = FALSE]
  }
  SC
}

# Pattern-wise casewise scores of the saturated (h1) Gaussian log-likelihood
# w.r.t. the sample statistics (mu, vech(Sigma)) under MISSING data (FIML).
# Rows = cases (all rows of Y, including those with NAs), columns = (means then
# vech(Sigma)) in psychonetrics ordering (the full p + p(p+1)/2 set; missing
# entries contribute exactly zero). For a case i with observed index set o and
# residual r = y_io - mu_o, W = Sigma_oo^-1 r:
#   mean part   : padded W   (zeros for missing variables)
#   cov  part   : 0.5 (W_a W_b - (Sigma_oo^-1)_ab) for a,b BOTH observed, else 0,
#                 with the vech (duplication) factor applied by doubling the
#                 off-diagonal entries (and leaving the diagonal as 0.5 * .).
# Cases are grouped by missing pattern for speed. Matches lavScores() under
# missing = "fiml" to machine precision (see extra fable_audit_scripts proto).
# FIML always carries a mean structure and is never a correlation input, so the
# meanstructure / corinput row dropping of ml_casewise_scores_h1 does not apply.
ml_casewise_scores_h1_missing <- function(Y, mu, Sigma){
  Y <- as.matrix(Y)
  n <- nrow(Y); p <- ncol(Y)
  ii <- unlist(lapply(seq_len(p), function(j) j:p))   # vech row idx (col-major)
  jj <- rep(seq_len(p), times = p:1)                  # vech col idx
  ns <- p * (p + 1) / 2
  SC <- matrix(0, n, p + ns)
  # Group cases by missing pattern (TRUE = observed):
  obsmat <- !is.na(Y)
  pat <- apply(obsmat, 1, function(z) paste(as.integer(z), collapse = ""))
  for (pt in unique(pat)){
    idx <- which(pat == pt)
    o <- which(obsmat[idx[1], ])
    if (length(o) == 0) next
    Soi <- solve(Sigma[o, o, drop = FALSE])
    Yc <- sweep(Y[idx, o, drop = FALSE], 2, mu[o])
    W <- Yc %*% Soi                       # rows: Sigma_oo^-1 (y_o - mu_o)
    SC[idx, o] <- W
    for (k in seq_len(ns)){
      a <- ii[k]; b <- jj[k]
      ia <- match(a, o); ib <- match(b, o)
      if (is.na(ia) || is.na(ib)) next    # contributes 0 unless both observed
      m <- 0.5 * (W[, ia] * W[, ib] - Soi[ia, ib])
      SC[idx, p + k] <- if (a == b) m else 2 * m
    }
  }
  SC
}

# Pattern-based OBSERVED information of the UNSTRUCTURED (saturated, h1) Gaussian
# log-likelihood w.r.t. (mu, vech(Sigma)) under MISSING data, evaluated at the
# supplied (mu, Sigma). This is the negative Hessian of the FIML saturated
# log-likelihood (lavaan's lav_mvnorm_missing_logl_hessian), divided by n so it
# is a UNIT information matrix. Used in the Yuan-Bentler-Mplus tr_h1 term under
# FIML. Per pattern p (frequency n_p, observed set o, pattern mean M_p, pattern
# covariance S_p with divisor n_p):
#   s_inv = pad(Sigma_oo^-1)
#   t21   = pad(Sigma_oo^-1 (M_p - mu_o))
#   Wt    = S_p + (M_p - mu_o)(M_p - mu_o)'
#   t22   = pad(Sigma_oo^-1 (2 Wt - Sigma_oo) Sigma_oo^-1)
# blocks [s_inv, t21-cross via D'(t21 (x) s_inv); 0.5 D'(s_inv (x) t22) D],
# frequency-weighted and divided by n. Matches lavaan's unstructured observed h1
# information under missingness to machine precision.
ml_h1_information_observed_missing <- function(Y, mu, Sigma){
  Y <- as.matrix(Y)
  n <- nrow(Y); p <- ncol(Y); ns <- p * (p + 1) / 2
  H11 <- matrix(0, p, p); H21 <- matrix(0, ns, p); H22 <- matrix(0, ns, ns)
  obsmat <- !is.na(Y)
  pat <- apply(obsmat, 1, function(z) paste(as.integer(z), collapse = ""))
  for (pt in unique(pat)){
    idx <- which(pat == pt); npat <- length(idx)
    o <- which(obsmat[idx[1], ])
    if (length(o) == 0) next
    Yp <- Y[idx, o, drop = FALSE]
    Mp <- colMeans(Yp)
    Sp <- crossprod(sweep(Yp, 2, Mp)) / npat
    Soi <- solve(Sigma[o, o, drop = FALSE])
    s_inv <- matrix(0, p, p); s_inv[o, o] <- Soi
    t21 <- matrix(0, p, 1); t21[o, 1] <- Soi %*% (Mp - mu[o])
    Wt <- Sp + tcrossprod(Mp - mu[o])
    aaa <- Soi %*% (2 * Wt - Sigma[o, o, drop = FALSE]) %*% Soi
    t22 <- matrix(0, p, p); t22[o, o] <- aaa
    H11 <- H11 + npat * s_inv
    H21 <- H21 + npat * lavaan::lav_matrix_duplication_pre(t21 %x% s_inv)
    H22 <- H22 + npat * 0.5 * lavaan::lav_matrix_duplication_pre_post(s_inv %x% t22)
  }
  rbind(cbind(H11, t(H21)), cbind(H21, H22)) / n
}

# Per-group EM/saturated (h1, unstructured) moments under FIML, as a list of
# list(mu, sigma) per group. These are a property of the DATA (not of the
# structural model) and match lavaan's lavTech(fit, "h1") moments. Source order:
#   1. the model's own fitted saturated submodel (x@baseline_saturated$saturated)
#      — available for a fitted structural model;
#   2. otherwise (e.g. the baseline model, which carries no saturated submodel),
#      fit a saturated varcov FIML model to the stored raw data per group.
# Returns NULL if neither source is available.
fiml_em_saturated_moments <- function(x){
  nGroups <- nrow(x@sample@groups)

  # 1. Reuse the model's fitted saturated submodel if present:
  sat <- x@baseline_saturated$saturated
  if (!is.null(sat) && is(sat, "psychonetrics") && length(sat@modelmatrices) >= nGroups){
    return(lapply(seq_len(nGroups), function(g)
      list(mu = as.numeric(sat@modelmatrices[[g]]$mu),
           sigma = as.matrix(sat@modelmatrices[[g]]$sigma))))
  }

  # 2. Fall back to fitting a saturated varcov FIML model to the raw data. The
  # saturated (fully free) varcov has df = 0, so its fitted moments are the EM
  # estimates of (mu, Sigma) under missingness.
  rawdata <- x@sample@rawdata
  vars <- attr(rawdata, "vars")
  groupcol <- attr(rawdata, "groups")
  if (is.null(vars) || is.null(groupcol) || nrow(rawdata) == 0) return(NULL)

  vc_args <- list(data = rawdata, vars = vars, estimator = "FIML",
                  storedata = FALSE, baseline_saturated = FALSE)
  if (nGroups > 1) vc_args$groups <- groupcol   # omit for a single group
  satmod <- tryCatch(do.call(varcov, vc_args), error = function(e) NULL)
  if (is.null(satmod)) return(NULL)
  satfit <- tryCatch(
    runmodel(satmod,
      addfit = FALSE, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE,
      analyticFisher = FALSE, verbose = FALSE, warn_gradient = FALSE,
      warn_bounds = FALSE, warn_improper = FALSE),
    error = function(e) NULL)
  if (is.null(satfit) || length(satfit@modelmatrices) < nGroups) return(NULL)

  # varcov() may reorder groups; align by the group labels used by x:
  sat_labels <- satfit@sample@groups$label
  x_labels <- x@sample@groups$label
  lapply(seq_len(nGroups), function(g){
    sg <- match(x_labels[g], sat_labels)
    if (is.na(sg)) sg <- g
    list(mu = as.numeric(satfit@modelmatrices[[sg]]$mu),
         sigma = as.matrix(satfit@modelmatrices[[sg]]$sigma))
  })
}

# Assemble the Huber-White (MLR) meat components for estimator = "ML" (complete
# data, Phase 1) or estimator = "FIML" (within-row missing data, Phase 2) from
# raw data. Returns NULL when raw data are unavailable.
# Components per group g:
#   Delta_g : group model Jacobian (means then vech(Sigma))
#   B1_g    : crossprod of casewise saturated scores at the MODEL-implied
#             moments (muHat_g, SigmaHat_g) / n_g  (the SE meat: casewise scores
#             w.r.t. theta are SC1_g %*% Delta_g). Under FIML the scores are the
#             pattern-wise FIML scores and n_g counts ALL rows (incl. incomplete).
#   B1u_g   : crossprod of casewise saturated scores at the SATURATED moments
#             / n_g  (for the Yuan-Bentler tr_h1 term). Complete data: the sample
#             moments (ybar_g, S_g). FIML: the EM/saturated moments (mu1_g,
#             Sigma1_g) from the saturated model's modelmatrices (= lavaan's h1).
#   S_g     : the moments used for the unstructured h1 information. Complete data:
#             sample covariance (divisor n_g, ML/biased — matches Gamma & lavaan).
#             FIML: the EM/saturated covariance Sigma1_g.
#   mu1_g, Y_g : (FIML only) the EM/saturated mean and the raw data matrix (with
#             NAs) for the group, needed by the pattern-based h1 observed
#             information in the scaled test.
#   missing : TRUE for the FIML (missing-data) path, FALSE for complete data.
#   fg      : group weights n_g / N
#   A       : full OBSERVED unit information of the model parameters
#             = 0.5 * jacobian(psychonetrics_gradient, parVector(x), x) symmetrised
#             (do NOT pre-apply expectedmodel(): that yields the EXPECTED info).
#             This is the observed information of the ML (complete) or FIML
#             (missing) fit function, which forms the Huber-White sandwich bread.
mlr_meat_components <- function(x){
  is_fiml <- x@estimator == "FIML"
  if (!is_fiml && x@estimator != "ML") return(NULL)
  if (length(x@modelmatrices) == 0) return(NULL)
  if (!.hasSlot(x@sample, "rawdata") || nrow(x@sample@rawdata) == 0) return(NULL)
  if (!is_fiml){
    if (length(x@sample@fimldata) > 0) return(NULL)   # complete-data path only
    if (!has_WLS_Gamma(x)) return(NULL)
  } else {
    if (length(x@sample@fimldata) == 0) return(NULL)  # FIML path needs patterns
  }

  rawdata <- x@sample@rawdata
  vars <- attr(rawdata, "vars")
  groupcol <- attr(rawdata, "groups")
  if (is.null(vars) || is.null(groupcol)) return(NULL)

  meanstructure <- x@meanstructure
  corinput <- isTRUE(x@sample@corinput)

  nGroups <- nrow(x@sample@groups)
  group_labels <- x@sample@groups$label
  nobs_per_group <- x@sample@groups$nobs
  nTotal <- sum(nobs_per_group)

  Delta_full <- build_Delta_full(x)

  # Row count per group (same ordering as Delta). For complete-data ML the
  # stored Gamma gives the per-group statistic count (which honours
  # meanstructure / corinput dropping); for FIML there is no Gamma and the
  # statistics are always the full means then vech(Sigma) per group.
  if (!is_fiml){
    nrows_per_group <- sapply(x@sample@WLS.Gamma, nrow)
  } else {
    nvar <- length(vars)
    nrows_per_group <- rep(nvar + nvar * (nvar + 1) / 2, nGroups)
  }
  if (sum(nrows_per_group) != nrow(Delta_full)) return(NULL)

  # Saturated (EM) moments per group for the FIML tr_h1 term (= lavaan's h1
  # unstructured moments). Reuses the model's saturated submodel when present,
  # otherwise fits a saturated varcov FIML to the raw data (e.g. for the baseline
  # model, which carries no saturated submodel of its own):
  sat_mm <- NULL
  if (is_fiml){
    sat_mm <- fiml_em_saturated_moments(x)
    if (is.null(sat_mm)) return(NULL)
  }

  Delta_list <- vector("list", nGroups)
  B1_list <- vector("list", nGroups)
  B1u_list <- vector("list", nGroups)
  S_list <- vector("list", nGroups)
  mu1_list <- if (is_fiml) vector("list", nGroups) else NULL
  Y_list <- if (is_fiml) vector("list", nGroups) else NULL

  row_offset <- 0
  for (g in seq_len(nGroups)){
    ng_stat <- nrows_per_group[g]
    row_inds <- row_offset + seq_len(ng_stat)
    Delta_list[[g]] <- Delta_full[row_inds, , drop = FALSE]
    row_offset <- row_offset + ng_stat

    Yg <- as.matrix(rawdata[rawdata[[groupcol]] == group_labels[g], vars, drop = FALSE])

    muHat <- as.numeric(x@modelmatrices[[g]]$mu)
    SigmaHat <- as.matrix(x@modelmatrices[[g]]$sigma)

    if (!is_fiml){
      # Complete data (Phase 1): drop incomplete rows (there are none here, but
      # keep the listwise guard for robustness) and use the normal-theory scores.
      Yg <- Yg[stats::complete.cases(Yg), , drop = FALSE]
      ng <- nrow(Yg)
      if (ng < 2) return(NULL)

      SC1 <- ml_casewise_scores_h1(Yg, muHat, SigmaHat,
                                   meanstructure = meanstructure, corinput = corinput)
      B1_list[[g]] <- crossprod(SC1) / ng

      ybar <- colMeans(Yg)
      Sg <- crossprod(sweep(Yg, 2, ybar)) / ng
      S_list[[g]] <- Sg
      SC1u <- ml_casewise_scores_h1(Yg, ybar, Sg,
                                    meanstructure = meanstructure, corinput = corinput)
      B1u_list[[g]] <- crossprod(SC1u) / ng
    } else {
      # FIML (Phase 2): keep ALL rows (incl. incomplete); pattern-wise scores.
      ng <- nrow(Yg)
      if (ng < 2) return(NULL)

      SC1 <- ml_casewise_scores_h1_missing(Yg, muHat, SigmaHat)
      B1_list[[g]] <- crossprod(SC1) / ng

      # Saturated (EM) moments for the unstructured h1 first-order info:
      mu1 <- as.numeric(sat_mm[[g]]$mu)
      Sigma1 <- as.matrix(sat_mm[[g]]$sigma)
      mu1_list[[g]] <- mu1
      S_list[[g]] <- Sigma1
      Y_list[[g]] <- Yg
      SC1u <- ml_casewise_scores_h1_missing(Yg, mu1, Sigma1)
      B1u_list[[g]] <- crossprod(SC1u) / ng
    }
  }

  # Full OBSERVED unit information of the model parameters (hessian-based).
  # numeric_FisherInformation() applies expectedmodel() first (giving the
  # EXPECTED information); here we deliberately do NOT, to obtain the observed
  # information that forms the Huber-White sandwich bread.
  G <- function(par) as.numeric(psychonetrics_gradient(par, x))
  A <- 0.5 * numDeriv::jacobian(G, parVector(x))
  A <- 0.5 * (A + t(A))

  list(
    Delta = Delta_list,
    B1 = B1_list,
    B1u = B1u_list,
    S = S_list,
    mu1 = mu1_list,    # FIML only (NULL otherwise): EM saturated means per group
    Y = Y_list,        # FIML only (NULL otherwise): raw data (with NAs) per group
    missing = is_fiml,
    A = A,
    fg = nobs_per_group / nTotal,
    nGroups = nGroups,
    npar = ncol(Delta_full),
    nTotal = nTotal,
    meanstructure = meanstructure,
    corinput = corinput
  )
}

# Assemble the per-group model Jacobian (Delta), weight matrix (W) and
# asymptotic covariance of the sample statistics (Gamma) for the least-squares
# estimators (WLS / DWLS / ULS). Both the W and Gamma lists follow the lavaan
# group-weighting convention: W is multiplied by f_g = n_g / n and Gamma is
# divided by f_g. With this scaling, sum_g Delta_g' W_g Delta_g equals the
# (unit-information) Fisher information used elsewhere, so the naive VCOV
# (1/n)(Delta'WDelta)^-1 and the robust sandwich
# (1/n)(Delta'WDelta)^-1 (Delta'W Gamma W Delta)(Delta'WDelta)^-1 are computed
# in consistent units (and coincide when W = Gamma^-1).
#
# Returns NULL when Gamma is unavailable (e.g. a model fit to a covariance
# matrix rather than raw data) or when the estimator is not least-squares.
# Reused by compute_wlsmv_correction() (scaled test statistic) and by
# getVCOV() (robust standard errors).
wls_sandwich_components <- function(x) {

  # Only defined for the least-squares estimators:
  if (!x@estimator %in% c("WLS", "DWLS", "ULS")) {
    return(NULL)
  }

  # Check that Gamma is available:
  if (!has_WLS_Gamma(x)) {
    return(NULL)
  }

  # Prepare model:
  if (x@cpp) {
    prep <- prepareModel_cpp(parVector(x), x)
  } else {
    prep <- prepareModel(parVector(x), x)
  }

  # Get model Jacobian (block-diagonal across groups):
  if (x@cpp) {
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov_cpp,
      "lvm" = d_phi_theta_lvm_cpp,
      "var1" = d_phi_theta_var1_cpp,
      "dlvm1" = d_phi_theta_dlvm1_cpp,
      "tsdlvm1" = d_phi_theta_tsdlvm1_cpp,
      "meta_varcov" = d_phi_theta_meta_varcov_cpp,
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising_cpp,
      "ml_lvm" = d_phi_theta_ml_lvm_cpp,
      "meta_lvm" = d_phi_theta_meta_lvm_cpp,
      "meta_var1" = d_phi_theta_meta_var1_cpp
    )
  } else {
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov,
      "lvm" = d_phi_theta_lvm,
      "var1" = d_phi_theta_var1,
      "dlvm1" = d_phi_theta_dlvm1,
      "tsdlvm1" = d_phi_theta_tsdlvm1,
      "meta_varcov" = d_phi_theta_meta_varcov,
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising,
      "ml_lvm" = d_phi_theta_ml_lvm,
      "meta_lvm" = d_phi_theta_meta_lvm,
      "meta_var1" = d_phi_theta_meta_var1
    )
  }

  modelPart <- modelJacobian(prep)

  # Get M matrix (maps constrained to full parameters):
  if (x@cpp) {
    M <- Mmatrix_cpp(x@parameters)
  } else {
    M <- Mmatrix(x@parameters)
  }

  # Full Delta = Jacobian %*% M (dsigma/dtheta_free)
  Delta_full <- as.matrix(modelPart %*% M)

  # Group info:
  nGroups <- nrow(x@sample@groups)
  nobs_per_group <- x@sample@groups$nobs
  nTotal <- sum(nobs_per_group)

  # Determine row ranges per group from WLS.W dimensions:
  nrows_per_group <- sapply(x@sample@WLS.W, nrow)

  # Build per-group W (scaled) and Gamma (scaled) lists, and extract Delta blocks:
  W_list <- vector("list", nGroups)
  Gamma_list <- vector("list", nGroups)
  Delta_list <- vector("list", nGroups)

  row_offset <- 0
  for (g in seq_len(nGroups)) {
    ng <- nrows_per_group[g]
    fg <- nobs_per_group[g] / nTotal  # group weight

    # Extract per-group Delta block:
    row_inds <- row_offset + seq_len(ng)
    Delta_list[[g]] <- Delta_full[row_inds, , drop = FALSE]
    row_offset <- row_offset + ng

    # Get W for this group (weight matrix used in estimation):
    W_g <- as.matrix(x@sample@WLS.W[[g]])
    if (x@estimator == "DWLS") {
      W_g <- diag(diag(W_g))
    } else if (x@estimator == "ULS") {
      W_g <- diag(nrow(W_g))
    }

    # Get Gamma for this group (full asymptotic covariance):
    Gamma_g <- as.matrix(x@sample@WLS.Gamma[[g]])

    # Scale by group weight (following lavaan convention):
    W_list[[g]] <- W_g * fg
    Gamma_list[[g]] <- Gamma_g / fg
  }

  # Bread: DtWD = sum_g Delta_g' W_g Delta_g (group-weighted). This equals the
  # (unit-information) Fisher information stored in x@information.
  npar_free <- ncol(Delta_full)
  DtWD <- matrix(0, npar_free, npar_free)
  for (g in seq_len(nGroups)) {
    DtWD <- DtWD + t(Delta_list[[g]]) %*% W_list[[g]] %*% Delta_list[[g]]
  }

  list(
    Delta = Delta_list,
    W = W_list,
    Gamma = Gamma_list,
    DtWD = DtWD,
    nGroups = nGroups,
    npar = npar_free,
    nTotal = nTotal
  )
}

# Integer degrees of freedom of a fitted model (nobs - npar, adjusted for a
# constrained saturated model). Shared by the scaled-test helpers.
model_df_integer <- function(x){
  df_integer <- x@sample@nobs - max(x@parameters$par)
  if (!is.null(x@baseline_saturated$saturated)) {
    df_integer <- df_integer -
      (x@baseline_saturated$saturated@sample@nobs -
         max(x@baseline_saturated$saturated@parameters$par))
  }
  df_integer
}

# Unscaled likelihood-ratio chi-square of a fitted ML model = -2*(LL - satLL).
# Used by the robust-ML scaled-test helpers when the caller (addfit) does not
# pass the already-computed chi-square. Returns NULL if a LL is unavailable.
# NB: for ML this is NOT objective * N (objective is the ML discrepancy).
ml_model_chisq <- function(x){
  LL <- tryCatch(psychonetrics_logLikelihood(x), error = function(e) NA_real_)
  satLL <- NA_real_
  if (!is.null(x@baseline_saturated$saturated)){
    satLL <- tryCatch(psychonetrics_logLikelihood(x@baseline_saturated$saturated),
                      error = function(e) NA_real_)
  }
  if (is.na(LL) || is.na(satLL)) return(NULL)
  -2 * (LL - satLL)
}

# Core Satorra-Bentler-family scaled test statistic computation, shared by the
# least-squares (WLS/DWLS/ULS) and robust-ML (MLM/MLMV/MLMVS) paths.
#
# Inputs (per group g, lavaan group-weighting convention already applied):
#   Delta_list  : list of group model Jacobians (means then vech(Sigma))
#   Wt_list     : list of WEIGHT blocks   Wtilde_g = f_g * W_g
#                 (W_g = the estimating weight matrix: V_g for ML, W_g for WLS,
#                  diag(W_g) for DWLS, I for ULS)
#   Gt_list     : list of GAMMA blocks    Gtilde_g = Gamma_g / f_g
#   bread       : the matrix E whose inverse projects out the fitted directions
#                 (DtWD = sum_g Delta_g' Wtilde_g Delta_g for least squares;
#                  the expected Fisher information x@information for ML — these
#                  coincide for ML since E = sum_g f_g Delta_g' V_g Delta_g)
#   chisq_naive : the (unscaled) chi-square = N * objective (or -2(LL-satLL))
#   df_integer  : integer degrees of freedom
#
# The U Gamma trace and especially trace((U Gamma)^2) are accumulated on the
# GLOBAL block-diagonal stack (V_all, G_all, Delta_all) rather than summed per
# group. The per-group accumulation is only correct when there are NO equality
# constraints linking parameters ACROSS groups; with cross-group constraints
# (e.g. group.equal in a multi-group model) U Gamma is not block-diagonal and
# the per-group sum of trace((U_g Gamma_g)^2) understates the true
# trace((U Gamma)^2). lavaan fixed the analogous bug in 0.6-13. For a single
# group the global form reduces exactly to the per-group form, so single-group
# output is unchanged.
compute_scaled_test_core <- function(Delta_list, Wt_list, Gt_list, bread,
                                      chisq_naive, df_integer){
  nGroups <- length(Delta_list)

  E_inv <- tryCatch(solve(bread), error = function(e) NULL)
  if (is.null(E_inv) || any(!is.finite(E_inv))) return(NULL)

  bdiag_dense <- function(L){
    n <- sum(vapply(L, nrow, integer(1)))
    M <- matrix(0, n, n)
    o <- 0L
    for (B in L){ nb <- nrow(B); M[o + seq_len(nb), o + seq_len(nb)] <- B; o <- o + nb }
    M
  }

  if (nGroups == 1L){
    # Single group: per-group form (identical to the global form, but cheaper).
    WD <- Wt_list[[1]] %*% Delta_list[[1]]
    U <- Wt_list[[1]] - WD %*% E_inv %*% t(WD)
    UG <- U %*% Gt_list[[1]]
  } else {
    V_all <- bdiag_dense(Wt_list)
    G_all <- bdiag_dense(Gt_list)
    Delta_all <- do.call(rbind, Delta_list)
    U <- V_all - V_all %*% Delta_all %*% E_inv %*% t(Delta_all) %*% V_all
    UG <- U %*% G_all
  }

  trace_UG  <- sum(diag(UG))
  trace_UG2 <- sum(UG * t(UG))   # tr((U Gamma)^2)

  if (trace_UG2 <= 0 || !is.finite(trace_UG2)) return(NULL)

  # Mean-and-Variance adjusted (Satterthwaite; fractional df):
  df_mv <- trace_UG^2 / trace_UG2
  scaling_mv <- trace_UG / df_mv  # = trace_UG2 / trace_UG
  chisq_mv <- chisq_naive / scaling_mv

  # Satorra-Bentler (mean adjusted): c = trace_UG / df; T = T_naive / c.
  scaling_sb <- trace_UG / df_integer
  chisq_sb <- chisq_naive / scaling_sb

  # Scaled-shifted (Satorra 2000; Asparouhov & Muthen 2010):
  a <- sqrt(df_integer / trace_UG2)
  shift <- df_integer - a * trace_UG
  chisq_ss <- a * chisq_naive + shift
  scaling_ss <- 1 / a

  list(
    # Satorra-Bentler (mean adjusted) — used by MLM:
    chisq.scaled.sb = chisq_sb,
    df.scaled.sb = df_integer,
    scaling.factor.sb = scaling_sb,
    # Mean-variance adjusted (fractional df) — used by MLMVS:
    chisq.scaled.mv = chisq_mv,
    df.scaled.mv = df_mv,
    scaling.factor.mv = scaling_mv,
    # Scaled-shifted (integer df, lavaan WLSMV / MLMV):
    chisq.scaled = chisq_ss,
    df.scaled = df_integer,
    scaling.factor = scaling_ss,
    shift.parameter = shift,
    # Raw traces for diagnostics:
    trace.UGamma = trace_UG,
    trace.UGamma2 = trace_UG2
  )
}

# Compute WLSMV (mean-and-variance / scaled-shifted) scaled test statistic for
# the least-squares estimators (WLS/DWLS/ULS).
# References: Satorra & Bentler (1994), Muthen (1993), Asparouhov & Muthen (2010)
compute_wlsmv_correction <- function(x) {

  # Assemble per-group Delta, W (scaled by f_g) and Gamma (scaled by 1/f_g):
  comp <- wls_sandwich_components(x)
  if (is.null(comp)) {
    return(NULL)
  }

  compute_scaled_test_core(
    Delta_list = comp$Delta,
    Wt_list = comp$W,        # already f_g * W_g
    Gt_list = comp$Gamma,    # already Gamma_g / f_g
    bread = comp$DtWD,       # sum_g Delta_g' (f_g W_g) Delta_g
    chisq_naive = x@objective * comp$nTotal,
    df_integer = model_df_integer(x)
  )
}

# Compute the robust-ML scaled test statistic (Satorra-Bentler family) for
# estimator = "ML" (used by MLM/MLMV/MLMVS). Uses the normal-theory weight
# matrix V_g and the expected Fisher information E = x@information.
# References: Satorra & Bentler (1994), Satorra (2000), Asparouhov & Muthen (2010)
compute_ml_scaled_test <- function(x, comp = NULL, chisq = NULL) {
  if (is.null(comp)) comp <- ml_robust_components(x)
  if (is.null(comp)) return(NULL)
  # The unscaled chi-square for ML is -2*(LL - satLL), NOT objective*N. The
  # caller (addfit) passes the chi-square it already computed; fall back to
  # computing it here if not supplied.
  if (is.null(chisq)) chisq <- ml_model_chisq(x)
  if (is.null(chisq) || !is.finite(chisq)) return(NULL)

  fg <- comp$fg
  nGroups <- comp$nGroups
  Wt_list <- vector("list", nGroups)
  Gt_list <- vector("list", nGroups)
  for (g in seq_len(nGroups)){
    Wt_list[[g]] <- fg[g] * comp$V[[g]]       # Vtilde_g = f_g V_g
    Gt_list[[g]] <- comp$Gamma[[g]] / fg[g]   # Gtilde_g = Gamma_g / f_g
  }

  # For ML the bread is the expected unit Fisher information E. By construction
  # E = sum_g f_g Delta_g' V_g Delta_g = sum_g Delta_g' Wtilde_g Delta_g, so it
  # plays exactly the role of DtWD in the least-squares path.
  compute_scaled_test_core(
    Delta_list = comp$Delta,
    Wt_list = Wt_list,
    Gt_list = Gt_list,
    bread = comp$E,
    chisq_naive = chisq,
    df_integer = model_df_integer(x)
  )
}

# Robust-ML Huber-White (MLR) Yuan-Bentler-Mplus scaled test statistic for
# estimator = "ML" (complete data, Phase 1) or "FIML" (missing data, Phase 2).
# References: Yuan & Bentler (2000); the "Mplus" variant matches lavaan's
#   test = "yuan.bentler.mplus" (default for estimator = "MLR").
#
# c = (tr_h1 - tr_h0) / df, with
#   tr_h1 = sum_g tr( B1u_g A1u_g^{-1} )   built at the SATURATED moments
#           (B1u_g = crossprod of casewise saturated scores at those moments
#            / n_g). Complete data: A1u_g = normal-theory expected info at the
#            sample moments (ybar_g, S_g). FIML: A1u_g = the pattern-based
#            UNSTRUCTURED observed h1 information at the EM moments (mu1_g,
#            Sigma1_g), which is the correct unstructured information under
#            missingness (the complete-data closed form does not apply).
#   tr_h0 = sum_g f_g tr( B0_g A^{-1} )    with B0_g STRUCTURED (the same meat
#           that feeds the Huber-White SEs) and A the full hessian-based
#           OBSERVED information of the model parameters.
compute_mlr_scaled_test <- function(x, meat = NULL, chisq = NULL) {
  if (!(x@estimator %in% c("ML", "FIML"))) return(NULL)
  if (is.null(meat)) meat <- mlr_meat_components(x)
  if (is.null(meat)) return(NULL)
  if (is.null(chisq)) chisq <- ml_model_chisq(x)
  if (is.null(chisq) || !is.finite(chisq)) return(NULL)

  df_integer <- model_df_integer(x)
  if (df_integer <= 0) return(NULL)

  is_fiml <- isTRUE(meat$missing)

  # tr_h1: unstructured (saturated) first-order vs h1 information at the
  # saturated moments, accumulated per group.
  tr_h1 <- 0
  for (g in seq_len(meat$nGroups)){
    if (!is_fiml){
      # Complete data: normal-theory expected info at the sample moments.
      A1u <- ml_Vmat(meat$S[[g]], meanstructure = meat$meanstructure,
                     corinput = meat$corinput)
    } else {
      # FIML: pattern-based observed h1 information at the EM moments.
      A1u <- ml_h1_information_observed_missing(meat$Y[[g]], meat$mu1[[g]],
                                                meat$S[[g]])
    }
    B1u <- meat$B1u[[g]]
    A1u_inv <- tryCatch(solve(A1u), error = function(e) NULL)
    if (is.null(A1u_inv)) return(NULL)
    tr_h1 <- tr_h1 + sum(B1u * t(A1u_inv))
  }

  # tr_h0: structured first-order information (the SE meat) against the full
  # observed information A.
  A_inv <- tryCatch(solve(meat$A), error = function(e) NULL)
  if (is.null(A_inv)) return(NULL)
  tr_h0 <- 0
  for (g in seq_len(meat$nGroups)){
    B0_g <- t(meat$Delta[[g]]) %*% meat$B1[[g]] %*% meat$Delta[[g]]
    tr_h0 <- tr_h0 + meat$fg[g] * sum(B0_g * t(A_inv))
  }

  trUG <- tr_h1 - tr_h0
  if (!is.finite(trUG) || trUG <= 0) return(NULL)
  scaling <- trUG / df_integer

  list(
    chisq.scaled = chisq / scaling,
    df.scaled = df_integer,
    scaling.factor = scaling,
    shift.parameter = NA_real_,
    trace.UGamma = trUG
  )
}

# ML discrepancy (Gaussian) of a model with implied moments (muHat, SigmaHat)
# against target moments (M, S), times n_g, summed over groups: the
# "complete-data" chi-square obtained by fitting the structural model (held at
# its FIML estimates) to the EM-completed sample moments. Used by the FIML-C
# robust fit indices (XX3 for the model, XX3.null for the baseline).
#   F_ML = log|Sigma| - log|S| + tr(S Sigma^-1) - p + (M - mu)' Sigma^-1 (M - mu)
fiml_completed_chisq <- function(implied_list, target_list, nobs_per_group){
  XX <- 0
  for (g in seq_along(implied_list)){
    muHat <- implied_list[[g]]$mu
    SigmaHat <- implied_list[[g]]$sigma
    M <- target_list[[g]]$mu
    S <- target_list[[g]]$sigma
    p <- ncol(SigmaHat)
    Si <- solve(SigmaHat)
    dmu <- M - muHat
    Fml <- as.numeric(determinant(SigmaHat, logarithm = TRUE)$modulus) -
      as.numeric(determinant(S, logarithm = TRUE)$modulus) +
      sum(S * t(Si)) - p + as.numeric(crossprod(dmu, Si %*% dmu))
    XX <- XX + nobs_per_group[g] * Fml
  }
  XX
}

# Robust FIML fit indices: FIML-Corrected RMSEA / CFI / TLI (Savalei 2010,
# "FIML-C(V3)"), matching lavaan's defaults for estimator = "MLR",
# missing = "fiml". The standard (complete-data) robust RMSEA/CFI/TLI formula
# does NOT hold under missingness; lavaan substitutes a corrected chi-square XX3
# (the model fit to the EM-completed moments), corrected df3 and a corrected
# scaling c.hat3 built from the missing-data information.
#
# Per group g, with observed set defined by the missingness patterns, EM moments
# (mu1_g, Sigma1_g), model-implied moments (muHat_g, SigmaHat_g), Jacobian
# Delta_g = d(mu, vech Sigma)/dtheta, f_g = n_g / n:
#   Wm_g  = pattern-based OBSERVED h1 information under missingness at (mu1, Sigma1)
#   Wc_g  = normal-theory (complete-data) information at the EM Sigma1
#   Jm_g  = first-order (crossprod of casewise FIML scores) at (mu1, Sigma1)
#   Gamma_g = Wm_g^{-1} Jm_g Wm_g^{-1}
# Stacked block-diagonal (lavaan weighting: Wc.g = f_g Wc, Wmi.g = f_g^-1 Wm^-1,
# Jm.g = f_g Jm, Gamma.f = f_g^-1 Gamma), with E.inv = (Delta' (f_g Wm) Delta)^-1
# the h1-information-based structured information inverse, the correction factor is
#   k = tr(Wc Gamma) - 2 tr(Delta' Jm Wm^-1 Wc Delta E.inv)
#         + tr(Delta' Jm Delta E.inv  (Delta' Wc Delta E.inv)')
#   c.hat3 = k / df3
# and XX3 = sum_g n_g F_ML(model implied vs EM moments). The baseline (XX3.null,
# c.hat3.null) uses the same Wc/Jm/Wm but the baseline Jacobian. Returns NULL if
# any building block is unavailable (so the caller omits the robust indices
# rather than emitting incorrect values).
# References: Savalei (2010); Zhang & Savalei (2022, FIML-C).
compute_mlr_fimlc_robust <- function(x, meat = NULL){
  if (x@estimator != "FIML") return(NULL)
  if (is.null(meat)) meat <- mlr_meat_components(x)
  if (is.null(meat) || !isTRUE(meat$missing)) return(NULL)

  df3 <- model_df_integer(x)
  if (df3 <= 0) return(NULL)

  nGroups <- meat$nGroups
  fg <- meat$fg
  nobs_per_group <- x@sample@groups$nobs

  # Per-group Wm / Wc / Jm and the model/EM implied moments:
  bdiag_dense <- function(L){
    n <- sum(vapply(L, nrow, integer(1)))
    M <- matrix(0, n, n); o <- 0L
    for (B in L){ nb <- nrow(B); M[o + seq_len(nb), o + seq_len(nb)] <- B; o <- o + nb }
    M
  }
  Wc.g <- Wmi.g <- Jm.g <- Wm.g <- Gamma.g <- vector("list", nGroups)
  implied_list <- em_list <- vector("list", nGroups)
  ok <- TRUE
  for (g in seq_len(nGroups)){
    mu1 <- meat$mu1[[g]]; Sigma1 <- meat$S[[g]]; Yg <- meat$Y[[g]]
    Wm <- tryCatch(ml_h1_information_observed_missing(Yg, mu1, Sigma1), error = function(e) NULL)
    Wc <- tryCatch(ml_Vmat(Sigma1, meanstructure = TRUE, corinput = FALSE), error = function(e) NULL)
    if (is.null(Wm) || is.null(Wc)){ ok <- FALSE; break }
    Wmi <- tryCatch(solve(Wm), error = function(e) NULL)
    if (is.null(Wmi)){ ok <- FALSE; break }
    Jm <- meat$B1u[[g]]    # crossprod of casewise FIML scores at (mu1, Sigma1) / n_g
    Wc.g[[g]] <- fg[g] * Wc
    Wm.g[[g]] <- fg[g] * Wm
    Wmi.g[[g]] <- (1 / fg[g]) * Wmi
    Jm.g[[g]] <- fg[g] * Jm
    Gamma.g[[g]] <- (1 / fg[g]) * (Wmi %*% Jm %*% Wmi)
    implied_list[[g]] <- list(mu = as.numeric(x@modelmatrices[[g]]$mu),
                              sigma = as.matrix(x@modelmatrices[[g]]$sigma))
    em_list[[g]] <- list(mu = mu1, sigma = Sigma1)
  }
  if (!ok) return(NULL)

  Wc.all <- bdiag_dense(Wc.g)
  Wmi.all <- bdiag_dense(Wmi.g)
  Jm.all <- bdiag_dense(Jm.g)
  Gamma.all <- bdiag_dense(Gamma.g)
  Wm.all <- bdiag_dense(Wm.g)
  Delta.all <- do.call(rbind, meat$Delta)

  E.inv <- tryCatch(solve(t(Delta.all) %*% Wm.all %*% Delta.all), error = function(e) NULL)
  if (is.null(E.inv)) return(NULL)

  fimlc_k <- function(Delta.all, E.inv){
    tr11 <- sum(Wc.all * Gamma.all)
    tr12 <- sum((t(Delta.all) %*% Jm.all %*% Wmi.all %*% Wc.all %*% Delta.all) * E.inv)
    DWcDE <- t(Delta.all) %*% Wc.all %*% Delta.all %*% E.inv
    tr22 <- sum((t(Delta.all) %*% Jm.all %*% Delta.all %*% E.inv) * t(DWcDE))
    tr11 - 2 * tr12 + tr22
  }

  k <- fimlc_k(Delta.all, E.inv)
  c.hat3 <- k / df3
  XX3 <- fiml_completed_chisq(implied_list, em_list, nobs_per_group)
  if (!is.finite(XX3) || !is.finite(c.hat3)) return(NULL)

  out <- list(XX3 = XX3, df3 = df3, c.hat3 = c.hat3,
              XX3.null = NA_real_, df3.null = NA_real_, c.hat3.null = NA_real_)

  # Baseline (independence) model: same Wc/Jm/Wm (data-level) but the baseline
  # Jacobian and implied moments. Needed for cfi.robust / tli.robust.
  bl <- x@baseline_saturated$baseline
  if (!is.null(bl) && is(bl, "psychonetrics") && bl@computed &&
      length(bl@modelmatrices) >= nGroups){
    DeltaB.all <- tryCatch(do.call(rbind, {
      DB <- build_Delta_full(bl)
      # split into the same per-group statistic blocks as the model:
      nstat <- vapply(meat$Delta, nrow, integer(1))
      off <- 0L
      lapply(seq_len(nGroups), function(g){
        rows <- off + seq_len(nstat[g]); off <<- off + nstat[g]
        DB[rows, , drop = FALSE]
      })
    }), error = function(e) NULL)
    if (!is.null(DeltaB.all)){
      E.invB <- tryCatch(solve(t(DeltaB.all) %*% Wm.all %*% DeltaB.all),
                          error = function(e) NULL)
      df3.null <- bl@sample@nobs - max(bl@parameters$par)
      if (!is.null(x@baseline_saturated$saturated)){
        df3.null <- df3.null - (x@baseline_saturated$saturated@sample@nobs -
                                  max(x@baseline_saturated$saturated@parameters$par))
      }
      bl_implied <- lapply(seq_len(nGroups), function(g)
        list(mu = as.numeric(bl@modelmatrices[[g]]$mu),
             sigma = as.matrix(bl@modelmatrices[[g]]$sigma)))
      if (!is.null(E.invB) && df3.null > 0){
        kb <- fimlc_k(DeltaB.all, E.invB)
        XX3.null <- fiml_completed_chisq(bl_implied, em_list, nobs_per_group)
        if (is.finite(XX3.null) && is.finite(kb)){
          out$XX3.null <- XX3.null
          out$df3.null <- df3.null
          out$c.hat3.null <- kb / df3.null
        }
      }
    }
  }

  out
}

# Dispatch the robust-ML scaled test statistic for a fitted model according to
# its robust configuration (set by setestimator / the constructors). Returns a
# normalised list with components chisq.scaled, df.scaled, scaling.factor,
# shift.parameter (NA when not scaled-shifted) and scaling.factor.sb (the
# mean-adjusted Satorra-Bentler factor, used by the robust fit indices for ALL
# MLM-family variants), or NULL if the test cannot be computed.
#   MLM   -> Satorra-Bentler (mean adjusted)
#   MLMV  -> scaled-and-shifted
#   MLMVS -> mean-and-variance adjusted (Satterthwaite, fractional df)
#   MLR   -> Yuan-Bentler-Mplus
robust_scaled_test <- function(x, chisq = NULL){
  cfg <- get_robust_config(x)
  test <- cfg$test
  if (is.null(test) || !nzchar(test)) return(NULL)

  if (test == "yuan.bentler.mplus"){
    res <- compute_mlr_scaled_test(x, chisq = chisq)
    if (is.null(res)) return(NULL)
    res$scaling.factor.sb <- res$scaling.factor   # YB scaling is already "mean adjusted"
    return(res)
  }

  # MLM family (robust.sem): compute all three scalings, then select.
  full <- compute_ml_scaled_test(x, chisq = chisq)
  if (is.null(full)) return(NULL)

  if (test == "satorra.bentler"){
    res <- list(
      chisq.scaled = full$chisq.scaled.sb,
      df.scaled = full$df.scaled.sb,
      scaling.factor = full$scaling.factor.sb,
      shift.parameter = NA_real_
    )
  } else if (test == "scaled.shifted"){
    res <- list(
      chisq.scaled = full$chisq.scaled,        # scaled-shifted
      df.scaled = full$df.scaled,
      scaling.factor = full$scaling.factor,    # 1/a
      shift.parameter = full$shift.parameter
    )
  } else if (test == "mean.var.adjusted"){
    res <- list(
      chisq.scaled = full$chisq.scaled.mv,
      df.scaled = full$df.scaled.mv,           # fractional df*
      scaling.factor = full$scaling.factor.mv,
      shift.parameter = NA_real_
    )
  } else {
    return(NULL)
  }
  # The robust RMSEA/CFI/TLI use the Satorra-Bentler (mean-adjusted) scaling
  # factor c regardless of which test variant is reported:
  res$scaling.factor.sb <- full$scaling.factor.sb
  res$trace.UGamma <- full$trace.UGamma
  res$trace.UGamma2 <- full$trace.UGamma2
  res
}


# Compute analytical saturated log-likelihood for FIML estimation.
# Avoids optimizing the (potentially huge) saturated varcov model by
# computing the pattern-specific ML directly from sample statistics.
# Formula per pattern: n_p/2 * [-log|S_p| - p_p*log(2*pi) - p_p]
fiml_saturated_loglikelihood <- function(x) {
  fimldata <- x@sample@fimldata
  nGroups <- nrow(x@sample@groups)

  total_ll <- 0

  for (g in seq_len(nGroups)) {
    patterns <- fimldata[[g]]

    for (p in seq_along(patterns)) {
      pat <- patterns[[p]]
      n_p <- pat$n
      S_p <- pat$S
      p_p <- sum(pat$obs)

      if (n_p < 1 || p_p < 1) next

      if (n_p < 2) {
        # Single observation: no covariance information, only constant
        ll_p <- n_p * (-p_p * log(2 * pi))
      } else {
        # Compute log pseudo-determinant (handles singular S_p)
        ev <- eigen(S_p, symmetric = TRUE, only.values = TRUE)$values
        pos_ev <- ev[ev > .Machine$double.eps * max(1, ev[1]) * p_p]

        if (length(pos_ev) == 0) {
          ll_p <- n_p * (-p_p * log(2 * pi))
        } else {
          log_det_S <- sum(log(pos_ev))
          rank_S <- length(pos_ev)
          ll_p <- n_p * (-log_det_S - p_p * log(2 * pi) - rank_S)
        }
      }

      total_ll <- total_ll + ll_p
    }
  }

  total_ll / 2
}

# Per-group Gaussian ML discrepancy F_g of a fitted varcov / lvm model: the
# normal-theory fit function evaluated at the model-implied moments against the
# stored sample moments,
#   F_g = log|Sigma| - log|S| + tr(S Sigma^-1) - p  [ + (ybar - mu)' Sigma^-1 (ybar - mu) ]
# (the mean term is included only when the model carries a mean structure). The
# sum_g n_g F_g equals -2*(LL - satLL); the Wishart chi-square replaces the n_g
# multipliers by (n_g - 1). Returns a numeric vector (one F_g per group), or NULL
# if the implied / sample moments are unavailable.
gaussian_discrepancy_per_group <- function(x){
  if (length(x@modelmatrices) == 0) return(NULL)
  if (length(x@sample@covs) == 0) return(NULL)
  nGroup <- nrow(x@sample@groups)
  meanstructure <- x@meanstructure
  Fvec <- numeric(nGroup)
  for (g in seq_len(nGroup)){
    Sigma <- as.matrix(x@modelmatrices[[g]]$sigma)
    S <- as.matrix(x@sample@covs[[g]])
    p <- ncol(Sigma)
    Si <- solve(Sigma)
    Fml <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) -
      as.numeric(determinant(S, logarithm = TRUE)$modulus) +
      sum(S * t(Si)) - p
    if (meanstructure){
      mu <- as.numeric(x@modelmatrices[[g]]$mu)
      ybar <- as.numeric(x@sample@means[[g]])
      if (length(mu) == p && length(ybar) == p && !anyNA(ybar)){
        dmu <- ybar - mu
        Fml <- Fml + as.numeric(crossprod(dmu, Si %*% dmu))
      }
    }
    Fvec[g] <- Fml
  }
  Fvec
}

# Wishart chi-square of a fitted Gaussian model = sum_g (n_g - 1) * F_g (vs the
# normal-theory chi-square sum_g n_g * F_g = -2*(LL - satLL)). Used for the model
# and the baseline under likelihood = "wishart". Returns NULL if F_g unavailable.
wishart_chisq <- function(x){
  Fvec <- gaussian_discrepancy_per_group(x)
  if (is.null(Fvec)) return(NULL)
  sum((x@sample@groups$nobs - 1) * Fvec)
}

# Computes fit measures
addfit <- function(
 x, #, ebicTuning = 0.25
 verbose
){
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  if (verbose){
    message("Adding fit measures...")
  }
  
  # If not computed, stop:
  if (!x@computed){
    stop("Model has not yet been run. Use runmodel(object) first!")
  }


  # Sample size:
  sampleSize <- sum(x@sample@groups$nobs)
  
  # Fitmeasures list:
  fitMeasures <- list()

  # Helper to add likelihood-based information criteria (independent of the
  # chisq/baseline machinery). Used both in the normal path and before the
  # early-exit returns below, so that logl/npar/df/AIC/BIC/ebic* are retained
  # even when chisq or the baseline model are not usable.
  add_ll_ic <- function(fm){
    LL <- fm$logl
    npar <- fm$npar
    nV <- fm$nvar
    fm$aic.ll  <- -2*LL + 2*npar
    fm$aic.ll2 <- -2*LL + 2*npar + (2*npar^2 + 2*npar)/(sampleSize - npar - 1)
    fm$bic  <- -2*LL + npar * log(sampleSize)
    N.star  <- (sampleSize + 2) / 24
    fm$bic2 <- -2*LL + npar * log(N.star)
    fm$ebic.25 <- -2*LL + npar * log(sampleSize) + 4 * npar * 0.25 * log(nV)
    fm$ebic.5  <- -2*LL + npar * log(sampleSize) + 4 * npar * 0.5  * log(nV)
    fm$ebic.75 <- -2*LL + npar * log(sampleSize) + 4 * npar * 0.75 * log(nV)
    fm$ebic1   <- -2*LL + npar * log(sampleSize) + 4 * npar * 1    * log(nV)
    fm
  }

  
  # log likelihoods:
  # Saturated:
  if (x@estimator %in% c("FIML","ML")){
    sat_method <- x@baseline_saturated$satMethod
    if (is.null(sat_method)) sat_method <- "default"

    if (sat_method == "analytic" && length(x@sample@fimldata) > 0) {
      # Forced analytical saturated LL:
      satLL <- fiml_saturated_loglikelihood(x)
    } else if (!is.null(x@baseline_saturated$saturated)){
      satLL <- psychonetrics_logLikelihood(x@baseline_saturated$saturated)
    } else {
      satLL <- NA
    }
    # Baseline:
    if (!is.null(x@baseline_saturated$baseline)){
      basLL <- psychonetrics_logLikelihood(x@baseline_saturated$baseline)
    } else {
      basLL <- NA
    }
    
    # Model:
    LL <-  psychonetrics_logLikelihood(x)

    # Fallback: if model-based saturated LL < model LL (optimizer failed),
    # use analytical saturated LL for FIML and warn:
    if (!is.na(satLL) && !is.na(LL) && satLL < LL &&
        sat_method != "analytic" &&
        x@estimator == "FIML" && length(x@sample@fimldata) > 0) {
      satLL_fallback <- fiml_saturated_loglikelihood(x)
      if (satLL_fallback >= LL) {
        warning("Saturated model optimization did not converge properly (saturated LL < model LL). ",
                "Using analytical saturated LL instead. Consider using saturated = 'analytic' in runmodel().")
        satLL <- satLL_fallback
      }
    }
  } else {
    satLL <- NA
    basLL <- NA
    LL <- NA
  }

  # fixed.x: the x-block is conditioned on rather than modelled, so the reported
  # log-likelihoods are the CONDITIONAL log-likelihoods (likelihood of the
  # endogenous variables given x), obtained by subtracting the saturated marginal
  # x-block log-likelihood from the joint log-likelihoods. This matches lavaan's
  # fixed.x reporting. The chi-square (-2*(LL - satLL)) is unaffected because the
  # same marginal cancels (the x-block is fixed identically in the model and the
  # saturated reference), so it is computed from the unadjusted joint LL/satLL
  # below; only the reported logl / AIC / BIC change.
  fx_marg <- 0
  if (length(get_fixed_x(x)$idx) > 0){
    fx_marg <- tryCatch(fixed_x_marginal_loglik(x), error = function(e) 0)
  }

  # Add to list:
  fitMeasures$logl <- LL - fx_marg
  fitMeasures$unrestricted.logl <- satLL - fx_marg
  fitMeasures$baseline.logl <- basLL - fx_marg

  # Number of variables:
  fitMeasures$nvar <- nVar <- nrow(x@sample@variables)
  
  # Number of observations:
  fitMeasures$nobs <- x@sample@nobs
    
  # Number of parameters:
  fitMeasures$npar <- max(x@parameters$par)
  
  # Degrees of freedom:
  fitMeasures$df <- fitMeasures$nobs - fitMeasures$npar    
  if (!is.null(x@baseline_saturated$saturated)){
    fitMeasures$df <- fitMeasures$df  - (x@baseline_saturated$saturated@sample@nobs - max(x@baseline_saturated$saturated@parameters$par))
  } 

  # Compute Fmin:
  fitMeasures$objective <- x@objective
  
  # Likelihood ratio:
  if (x@estimator %in% c("FIML","ML")){
    # Under the Wishart Gaussian likelihood the chi-square uses (n_g - 1)
    # multipliers (chisq = sum_g (n_g-1) F_g) rather than the normal-theory
    # -2*(LL - satLL) = sum_g n_g F_g (matches lavaan likelihood = "wishart").
    wchisq <- if (is_wishart(x)) tryCatch(wishart_chisq(x), error = function(e) NULL) else NULL
    if (!is.null(wchisq) && is.finite(wchisq)){
      fitMeasures$chisq <- wchisq
    } else {
      fitMeasures$chisq <- -2 * (LL - satLL)
    }
  } else  if (x@estimator %in% c("WLS","DWLS","ULS")){
    fitMeasures$chisq <- x@objective  * (sampleSize)
  }
  fitMeasures$pvalue <- pchisq(fitMeasures$chisq, fitMeasures$df, lower.tail = FALSE)

  # WLSMV scaled test statistic (mean-and-variance adjusted):
  wlsmv_model <- NULL
  if (x@estimator %in% c("WLS","DWLS","ULS") && has_WLS_Gamma(x)) {
    experimentalWarning("WLSMV scaled test statistic")
    wlsmv_model <- tryCatch(compute_wlsmv_correction(x), error = function(e) NULL)
    if (!is.null(wlsmv_model)) {
      fitMeasures$chisq.scaled <- wlsmv_model$chisq.scaled
      fitMeasures$df.scaled <- wlsmv_model$df.scaled
      fitMeasures$pvalue.scaled <- pchisq(wlsmv_model$chisq.scaled, wlsmv_model$df.scaled, lower.tail = FALSE)
      fitMeasures$chisq.scaling.factor <- wlsmv_model$scaling.factor
    }
  }

  # Robust ML scaled test statistic (MLM/MLMV/MLMVS = Satorra-Bentler family;
  # MLR = Yuan-Bentler-Mplus). Stored alongside the unscaled chisq. MLR with
  # within-row missing data maps internally to estimator = "FIML":
  robust_model <- NULL
  if (x@estimator %in% c("ML", "FIML") && is_robust_ML(x)) {
    robust_model <- tryCatch(robust_scaled_test(x, chisq = fitMeasures$chisq), error = function(e) NULL)
    if (!is.null(robust_model)) {
      fitMeasures$chisq.scaled <- robust_model$chisq.scaled
      fitMeasures$df.scaled <- robust_model$df.scaled
      fitMeasures$pvalue.scaled <- pchisq(robust_model$chisq.scaled, robust_model$df.scaled, lower.tail = FALSE)
      fitMeasures$chisq.scaling.factor <- robust_model$scaling.factor
      if (!is.null(robust_model$shift.parameter) && is.finite(robust_model$shift.parameter)){
        fitMeasures$chisq.shift.parameters <- robust_model$shift.parameter
      }
    }
  }

  # FIML-Corrected robust fit indices (Savalei 2010, FIML-C(V3)) for MLR with
  # missing data. The standard (complete-data) robust RMSEA/CFI/TLI formula does
  # not apply under missingness; this returns the corrected XX3/df3/c.hat3 (and
  # the baseline analogues) used by the robust RMSEA/CFI/TLI blocks below. NULL
  # for complete data or when a building block is unavailable, in which case the
  # standard robust-index path is used (and, under FIML, the robust indices are
  # simply omitted rather than computed incorrectly).
  fimlc_robust <- NULL
  if (!is.null(robust_model) && x@estimator == "FIML" && is_robust_ML(x) &&
      identical(get_robust_config(x)$label, "MLR")) {
    fimlc_robust <- tryCatch(compute_mlr_fimlc_robust(x), error = function(e) NULL)
  }

  # Some pars:
  Tm <- fitMeasures$chisq
  dfm <- fitMeasures$df
  
  # Baseline model:
  if (!is.null(x@baseline_saturated$baseline) && x@baseline_saturated$baseline@computed){
    if (length(x@baseline_saturated$baseline@objective) == 0){
      x@baseline_saturated$baseline@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$baseline),x@baseline_saturated$baseline)
    }
    
    # fitMeasures$fmin_baseline <- x@baseline_saturated$baseline@objective
    # fitMeasures$baseline.chisq <-  sampleSize * fitMeasures$fmin_baseline
    if (x@estimator%in% c("FIML","ML")){
      # Wishart: baseline chisq also uses (n_g - 1) multipliers.
      wbchisq <- if (is_wishart(x)) tryCatch(wishart_chisq(x@baseline_saturated$baseline), error = function(e) NULL) else NULL
      if (!is.null(wbchisq) && is.finite(wbchisq)){
        fitMeasures$baseline.chisq <- wbchisq
      } else {
        fitMeasures$baseline.chisq <-  -2 * (basLL - satLL)
      }
    } else  if (x@estimator %in% c("WLS","DWLS","ULS")){
      fitMeasures$baseline.chisq <- x@baseline_saturated$baseline@objective  * (sampleSize)
    }
    fitMeasures$baseline.npar <- max(x@baseline_saturated$baseline@parameters$par)
    fitMeasures$baseline.df <- fitMeasures$nobs - max(x@baseline_saturated$baseline@parameters$par)
    if (!is.null(x@baseline_saturated$saturated)){
      # Mirror the adjustment applied to fitMeasures$df above: the baseline
      # chisq (-2 * (basLL - satLL)) is referenced against the saturated
      # model, so its df is npar_sat - npar_baseline. When the saturated
      # model itself is constrained (e.g. multi-group with equal=) so
      # npar_sat < nobs, baseline.df must subtract (nobs_sat - npar_sat).
      # If the saturated is fully saturated this subtracts zero (legacy case).
      fitMeasures$baseline.df <- fitMeasures$baseline.df - (x@baseline_saturated$saturated@sample@nobs - max(x@baseline_saturated$saturated@parameters$par))
    }
    fitMeasures$baseline.pvalue <- pchisq(fitMeasures$baseline.chisq, fitMeasures$baseline.df, lower.tail = FALSE)

    # WLSMV scaled baseline:
    wlsmv_baseline <- NULL
    if (!is.null(wlsmv_model) && x@estimator %in% c("WLS","DWLS","ULS") && has_WLS_Gamma(x)) {
      wlsmv_baseline <- tryCatch(compute_wlsmv_correction(x@baseline_saturated$baseline), error = function(e) NULL)
      if (!is.null(wlsmv_baseline)) {
        fitMeasures$baseline.chisq.scaled <- wlsmv_baseline$chisq.scaled
        fitMeasures$baseline.df.scaled <- wlsmv_baseline$df.scaled
        fitMeasures$baseline.pvalue.scaled <- pchisq(wlsmv_baseline$chisq.scaled, wlsmv_baseline$df.scaled, lower.tail = FALSE)
        fitMeasures$baseline.chisq.scaling.factor <- wlsmv_baseline$scaling.factor
      }
    }

    # Robust ML scaled baseline (same correction as the model; needed for the
    # robust incremental fit indices and the robust RMSEA baseline factor c_B):
    robust_baseline <- NULL
    if (!is.null(robust_model) && x@estimator %in% c("ML", "FIML") && is_robust_ML(x)) {
      robust_baseline <- tryCatch(
        robust_scaled_test(x@baseline_saturated$baseline, chisq = fitMeasures$baseline.chisq),
        error = function(e) NULL)
      if (!is.null(robust_baseline)) {
        fitMeasures$baseline.chisq.scaled <- robust_baseline$chisq.scaled
        fitMeasures$baseline.df.scaled <- robust_baseline$df.scaled
        fitMeasures$baseline.pvalue.scaled <- pchisq(robust_baseline$chisq.scaled, robust_baseline$df.scaled, lower.tail = FALSE)
        fitMeasures$baseline.chisq.scaling.factor <- robust_baseline$scaling.factor
      }
    }

    # Incremental Fit Indices
    Tb <- fitMeasures$baseline.chisq
    
    dfb <- fitMeasures$baseline.df
    # 
    # t1 <- (X2 - df)*df.null
    # t2 <- (X2.null - df.null)*df
    # if(df > 0 && abs(t2) > 0) {
    #   indices["tli"] <- indices["nnfi"] <- 1 - t1/t2
    # } else {
    #   indices["tli"] <- indices["nnfi"] <- 1
    # }
    
    fitMeasures$nfi <- (Tb - Tm) / Tb

    # Stop here if baseline is not good:
    if (is.null(dfb) || !is.finite(dfb) || !is.finite(Tb)){
      fitMeasures <- add_ll_ic(fitMeasures)
      x@fitmeasures <- fitMeasures
      return(x)
    }
    
    if(dfb > 0 && Tb > 0) {
      t1 <- Tb - Tm
      t2 <- Tb
      fitMeasures$pnfi <- (dfm/dfb) * t1/t2
    } else {
      fitMeasures$pnfi <- as.numeric(NA)
    }
    
    fitMeasures$tli <-  (Tb/dfb - Tm/dfm) / (Tb/dfb - 1) 

    
    t1 <- (Tm - dfm)*dfb
    t2 <- (Tb - dfb)*dfm
    fitMeasures$nnfi <- ifelse(dfm > 0 & abs(t2) > 0, 1 - t1/t2, 1)
    

    
    rfi_val <- (Tb/dfb - Tm/dfm) / (Tb/dfb)
    fitMeasures$rfi <- ifelse(Tb/dfb <= 0, 1, rfi_val)
    ifi_denom <- Tb - dfm
    fitMeasures$ifi <- ifelse(ifi_denom <= 0, 1, (Tb - Tm) / ifi_denom)
    fitMeasures$rni <-  ((Tb- dfb) - (Tm - dfm)) / (Tb - dfb)
    fitMeasures$cfi <- ifelse(dfm > Tm, 1, 1 - (Tm - dfm)/(Tb - dfb))

    # Scaled incremental fit indices (WLSMV):
    if (!is.null(wlsmv_model) && !is.null(wlsmv_baseline)) {
      Tm_s <- fitMeasures$chisq.scaled
      dfm_s <- fitMeasures$df.scaled
      Tb_s <- fitMeasures$baseline.chisq.scaled
      dfb_s <- fitMeasures$baseline.df.scaled

      fitMeasures$nfi.scaled <- (Tb_s - Tm_s) / Tb_s
      fitMeasures$tli.scaled <- (Tb_s/dfb_s - Tm_s/dfm_s) / (Tb_s/dfb_s - 1)
      t1_s <- (Tm_s - dfm_s)*dfb_s
      t2_s <- (Tb_s - dfb_s)*dfm_s
      fitMeasures$nnfi.scaled <- ifelse(dfm_s > 0 & abs(t2_s) > 0, 1 - t1_s/t2_s, 1)
      rfi_val_s <- (Tb_s/dfb_s - Tm_s/dfm_s) / (Tb_s/dfb_s)
      fitMeasures$rfi.scaled <- ifelse(Tb_s/dfb_s <= 0, 1, rfi_val_s)
      ifi_denom_s <- Tb_s - dfm_s
      fitMeasures$ifi.scaled <- ifelse(ifi_denom_s <= 0, 1, (Tb_s - Tm_s) / ifi_denom_s)
      fitMeasures$rni.scaled <- ((Tb_s - dfb_s) - (Tm_s - dfm_s)) / (Tb_s - dfb_s)
      fitMeasures$cfi.scaled <- ifelse(dfm_s > Tm_s, 1, 1 - (Tm_s - dfm_s)/(Tb_s - dfb_s))
    }

    # Robust incremental fit indices (Brosseau-Liard & Savalei 2014) for robust
    # ML. These use the UNscaled chi-squares Tm/Tb with the Satorra-Bentler
    # (mean-adjusted) scaling factors c (model) and c_B (baseline):
    #   cfi.robust = 1 - max(Tm - c dfm, 0) / max(Tm - c dfm, Tb - c_B dfb, 0)
    #   tli.robust = 1 - (Tm - c dfm) dfb / ((Tb - c_B dfb) dfm)
    # Under FIML (MLR + missing data) the corrected FIML-C(V3) chi-squares and
    # scaling factors (XX3/c.hat3 for the model, XX3.null/c.hat3.null for the
    # baseline) replace Tm/c and Tb/c_B (Savalei 2010), matching lavaan.
    if (!is.null(fimlc_robust) && is.finite(fimlc_robust$c.hat3) &&
        is.finite(fimlc_robust$c.hat3.null)) {
      Tm_r <- fimlc_robust$XX3;      dfm_r <- fimlc_robust$df3;      c_m <- fimlc_robust$c.hat3
      Tb_r <- fimlc_robust$XX3.null; dfb_r <- fimlc_robust$df3.null; c_b <- fimlc_robust$c.hat3.null
      t1r <- max(c(Tm_r - c_m*dfm_r, 0))
      t2r <- max(c(Tm_r - c_m*dfm_r, Tb_r - c_b*dfb_r, 0))
      fitMeasures$cfi.robust <- if (t2r <= 0) 1 else 1 - t1r/t2r
      ttr1 <- (Tm_r - c_m*dfm_r)*dfb_r
      ttr2 <- (Tb_r - c_b*dfb_r)*dfm_r
      fitMeasures$tli.robust <- if (dfm_r > 0 && abs(ttr2) > 0) 1 - ttr1/ttr2 else 1
      fitMeasures$nnfi.robust <- fitMeasures$tli.robust
    } else if (is.null(fimlc_robust) && x@estimator == "FIML" && is_robust_ML(x)) {
      # FIML-MLR but the FIML-C correction was unavailable: omit the robust
      # incremental indices rather than emit the (inapplicable) complete-data form.
    } else if (!is.null(robust_model) && !is.null(robust_baseline)) {
      c_m <- robust_model$scaling.factor.sb
      c_b <- robust_baseline$scaling.factor.sb
      if (!is.null(c_m) && !is.null(c_b) && is.finite(c_m) && is.finite(c_b)) {
        t1r <- max(c(Tm - c_m*dfm, 0))
        t2r <- max(c(Tm - c_m*dfm, Tb - c_b*dfb, 0))
        fitMeasures$cfi.robust <- if (t2r <= 0) 1 else 1 - t1r/t2r
        ttr1 <- (Tm - c_m*dfm)*dfb
        ttr2 <- (Tb - c_b*dfb)*dfm
        fitMeasures$tli.robust <- if (dfm > 0 && abs(ttr2) > 0) 1 - ttr1/ttr2 else 1
        fitMeasures$nnfi.robust <- fitMeasures$tli.robust
      }
    }

    } else {
    warning("No baseline model found, cannot add incremental fit indices...")
      fitMeasures$fmin_baseline <- NA
      fitMeasures$baseline.chisq <- NA
      fitMeasures$baseline.df <- NA
      fitMeasures$baseline.pvalue <- NA
      
      # Incremental Fit Indices
      fitMeasures$nfi <- NA
      fitMeasures$tli <- NA
      fitMeasures$nnfi <- NA
      fitMeasures$rfi <- NA
      fitMeasures$ifi <- NA
      fitMeasures$rni <- NA
      fitMeasures$cfi <- NA
  }
  
  # If LLs are not good, break here:

  if (!is.finite(Tm) ){
    fitMeasures <- add_ll_ic(fitMeasures)
    x@fitmeasures <- fitMeasures
    return(x)
  }
  
  # RMSEA

  # fitMeasures$rmsea <- sqrt( max(Tm - dfm,0) / (sampleSize * dfm))
  fitMeasures$rmsea <-  sqrt( max( c((Tm/sampleSize)/dfm - 1/sampleSize, 0) ) )
  if (!is.finite(fitMeasures$rmsea)) fitMeasures$rmsea <- NA
  
  # FIXME: Multi-group correction from lavaan source code:
  nGroups <- nrow(x@sample@groups)
  fitMeasures$rmsea <-  fitMeasures$rmsea  * sqrt(nGroups)

  # Guard against astronomically large chi-square values (e.g., diverged fits):
  # the noncentral chi-square searches below cannot be evaluated reliably then
  # and would flood the console with pnchisq convergence warnings:
  chisq_too_large <- is.finite(Tm) && Tm > 1e10

  # Codes for rmsea confidence interval taken from lavaan:
  lower.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.95)
  }
  if(is.na(Tm) || is.na(dfm) || chisq_too_large) {
    fitMeasures$rmsea.ci.lower <- NA
  } else if(dfm < 1 || lower.lambda(0) < 0.0) {
    fitMeasures$rmsea.ci.lower <- 0
  } else {
    if (lower.lambda(0) * lower.lambda(Tm) > 0){
      lambda.l <- NA
    } else {
      lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=Tm)$root,
                      silent=TRUE)
    }

    # Guard against a failed uniroot (mirror the upper-bound branch):
    if(inherits(lambda.l, "try-error")) {lambda.l <- NA }

    fitMeasures$rmsea.ci.lower <- sqrt( lambda.l/(sampleSize*dfm) ) * sqrt(nGroups)
  }
  
  N.RMSEA <- max(sampleSize, Tm*4) 
  upper.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.05)
  }
  if(is.na(Tm) || is.na(dfm) || chisq_too_large) {
    fitMeasures$rmsea.ci.upper <- NA
  } else if(dfm < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
    fitMeasures$rmsea.ci.upper <- 0
  } else {
    
    if (upper.lambda(0) * upper.lambda(N.RMSEA) > 0){
      lambda.u <- NA
    } else {
      
      lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                      silent=TRUE)  
    }
    
    if(inherits(lambda.u, "try-error")) {lambda.u <- NA }
    
    fitMeasures$rmsea.ci.upper <- sqrt( lambda.u/(sampleSize*dfm) )  * sqrt(nGroups)
  }
  
  fitMeasures$rmsea.pvalue <- if (chisq_too_large) NA else
    1 - pchisq(Tm, df=dfm, ncp=(sampleSize*dfm*0.05^2/nGroups))

  # Scaled RMSEA (WLSMV):
  if (!is.null(wlsmv_model)) {
    Tm_s <- fitMeasures$chisq.scaled
    dfm_s <- fitMeasures$df.scaled
    fitMeasures$rmsea.scaled <- sqrt( max( c((Tm_s/sampleSize)/dfm_s - 1/sampleSize, 0) ) )
    if (!is.finite(fitMeasures$rmsea.scaled)) fitMeasures$rmsea.scaled <- NA
    fitMeasures$rmsea.scaled <- fitMeasures$rmsea.scaled * sqrt(nGroups)
  }

  # Robust RMSEA (Brosseau-Liard & Savalei 2014; Savalei 2010) for robust ML.
  # Uses the UNscaled Tm with the Satorra-Bentler (mean-adjusted) scaling factor
  # c = trace(U Gamma)/df:
  #   rmsea.robust = sqrt( max( (Tm/N)/dfm - c/N, 0 ) ) * sqrt(G)
  # The same mean-adjusted c is used for ALL MLM-family variants (MLM/MLMV/MLMVS)
  # and matches lavaan's rmsea.robust (which is identical across these variants).
  # The confidence-interval bounds invert pchisq(Tm_scaled, dfm, ncp) for the
  # non-centrality lambda and map it back with sqrt(c lambda / (N dfm)) * sqrt(G).
  # Under FIML (MLR + missing data) the corrected FIML-C(V3) statistic XX3, df3
  # and scaling c.hat3 replace Tm, dfm and c (Savalei 2010), and the CI inverts
  # the noncentral chi-square on XX3.scaled = XX3 / c.hat3.
  use_fimlc_rmsea <- !is.null(fimlc_robust) && is.finite(fimlc_robust$c.hat3) &&
    fimlc_robust$df3 > 0
  if (use_fimlc_rmsea) {
    Tm_r <- fimlc_robust$XX3; dfm_r <- fimlc_robust$df3; c_rmsea <- fimlc_robust$c.hat3
    rr <- sqrt( max( c((Tm_r/sampleSize)/dfm_r - c_rmsea/sampleSize, 0) ) )
    if (!is.finite(rr)) rr <- NA
    fitMeasures$rmsea.robust <- rr * sqrt(nGroups)

    # CI on the scaled corrected statistic XX3.scaled = XX3 / c.hat3:
    Tm_sc <- Tm_r / c_rmsea
    dfm_sc <- dfm_r
    if (!chisq_too_large && is.finite(Tm_sc) && is.finite(dfm_sc) && dfm_sc >= 1) {
      low.l <- function(lambda) (pchisq(Tm_sc, df = dfm_sc, ncp = lambda) - 0.95)
      if (low.l(0) < 0.0){
        fitMeasures$rmsea.ci.lower.robust <- 0
      } else {
        ll <- if (low.l(0) * low.l(Tm_sc) > 0) NA else
          try(uniroot(low.l, lower = 0, upper = Tm_sc)$root, silent = TRUE)
        if (inherits(ll, "try-error")) ll <- NA
        fitMeasures$rmsea.ci.lower.robust <- sqrt(c_rmsea * ll / (sampleSize * dfm_r)) * sqrt(nGroups)
      }
      N.RMSEA.r <- max(sampleSize, Tm_sc * 4)
      up.l <- function(lambda) (pchisq(Tm_sc, df = dfm_sc, ncp = lambda) - 0.05)
      if (up.l(N.RMSEA.r) > 0 || up.l(0) < 0){
        fitMeasures$rmsea.ci.upper.robust <- 0
      } else {
        lu <- if (up.l(0) * up.l(N.RMSEA.r) > 0) NA else
          try(uniroot(up.l, lower = 0, upper = N.RMSEA.r)$root, silent = TRUE)
        if (inherits(lu, "try-error")) lu <- NA
        fitMeasures$rmsea.ci.upper.robust <- sqrt(c_rmsea * lu / (sampleSize * dfm_r)) * sqrt(nGroups)
      }
    }
  } else if (is.null(fimlc_robust) && x@estimator == "FIML" && is_robust_ML(x)) {
    # FIML-MLR but the FIML-C correction was unavailable: omit rmsea.robust
    # rather than emit the (inapplicable) complete-data form.
  } else if (!is.null(robust_model) && dfm > 0) {
    c_rmsea <- robust_model$scaling.factor.sb
    if (!is.null(c_rmsea) && is.finite(c_rmsea)) {
      rr <- sqrt( max( c((Tm/sampleSize)/dfm - c_rmsea/sampleSize, 0) ) )
      if (!is.finite(rr)) rr <- NA
      fitMeasures$rmsea.robust <- rr * sqrt(nGroups)

      # Confidence interval (invert the noncentral chi-square on the SCALED Tm):
      Tm_sc <- fitMeasures$chisq.scaled
      dfm_sc <- fitMeasures$df.scaled
      if (!chisq_too_large && is.finite(Tm_sc) && is.finite(dfm_sc) && dfm_sc >= 1) {
        low.l <- function(lambda) (pchisq(Tm_sc, df = dfm_sc, ncp = lambda) - 0.95)
        if (low.l(0) < 0.0){
          fitMeasures$rmsea.ci.lower.robust <- 0
        } else {
          ll <- if (low.l(0) * low.l(Tm_sc) > 0) NA else
            try(uniroot(low.l, lower = 0, upper = Tm_sc)$root, silent = TRUE)
          if (inherits(ll, "try-error")) ll <- NA
          fitMeasures$rmsea.ci.lower.robust <- sqrt(c_rmsea * ll / (sampleSize * dfm)) * sqrt(nGroups)
        }
        N.RMSEA.r <- max(sampleSize, Tm_sc * 4)
        up.l <- function(lambda) (pchisq(Tm_sc, df = dfm_sc, ncp = lambda) - 0.05)
        if (up.l(N.RMSEA.r) > 0 || up.l(0) < 0){
          fitMeasures$rmsea.ci.upper.robust <- 0
        } else {
          lu <- if (up.l(0) * up.l(N.RMSEA.r) > 0) NA else
            try(uniroot(up.l, lower = 0, upper = N.RMSEA.r)$root, silent = TRUE)
          if (inherits(lu, "try-error")) lu <- NA
          fitMeasures$rmsea.ci.upper.robust <- sqrt(c_rmsea * lu / (sampleSize * dfm)) * sqrt(nGroups)
        }
      }
    }
  }

  # SRMR (Bentler, 1995):
  # Standardized Root Mean Square Residual.
  # Standardize both observed and implied covariances by observed SDs.
  # When a mean structure is present, include standardized mean residuals
  # (matching lavaan's default srmr_bentler behavior).
  # Under FIML, use the saturated model's implied moments as "observed"
  # (EM estimate, matching lavaan's h1 approach).
  has_sigma <- length(x@modelmatrices) > 0 &&
    !is.null(x@modelmatrices[[1]]$sigma) &&
    is.matrix(x@modelmatrices[[1]]$sigma)
  if (length(x@sample@covs) > 0 && has_sigma) {
    use_saturated <- x@estimator == "FIML" &&
      !is.null(x@baseline_saturated$saturated) &&
      x@baseline_saturated$saturated@computed
    has_means <- x@meanstructure
    srmr_groups <- numeric(nGroups)
    nobs_per_group <- x@sample@groups$nobs
    for (g in seq_len(nGroups)) {
      if (use_saturated) {
        sat_sig <- x@baseline_saturated$saturated@modelmatrices[[g]]$sigma
        if (!is.matrix(sat_sig)) { srmr_groups[g] <- NA_real_; next }
        S <- as.matrix(sat_sig)
      } else {
        S <- as.matrix(x@sample@covs[[g]])
      }
      Sigma <- as.matrix(x@modelmatrices[[g]]$sigma)
      obs_sd <- sqrt(diag(S))
      # Guard against zero/near-zero variances:
      obs_sd[obs_sd < sqrt(.Machine$double.eps)] <- NA
      scale_mat <- tcrossprod(obs_sd)
      E <- (S - Sigma) / scale_mat
      e_cov <- E[lower.tri(E, diag = TRUE)]
      e_cov <- e_cov[!is.na(e_cov)]
      # Include mean residuals if mean structure is present:
      if (has_means) {
        if (use_saturated) {
          obs_mu <- as.numeric(x@baseline_saturated$saturated@modelmatrices[[g]]$mu)
        } else {
          obs_mu <- as.numeric(x@sample@means[[g]])
        }
        imp_mu <- as.numeric(x@modelmatrices[[g]]$mu)
        e_mean <- (obs_mu - imp_mu) / obs_sd
        e_mean <- e_mean[!is.na(e_mean)]
        e_all <- c(e_mean, e_cov)
      } else {
        e_all <- e_cov
      }
      srmr_groups[g] <- if (length(e_all) > 0) sqrt(mean(e_all^2)) else NA_real_
    }
    fitMeasures$srmr <- sum(nobs_per_group / sampleSize * srmr_groups, na.rm = TRUE)
    if (all(is.na(srmr_groups))) fitMeasures$srmr <- NA
  } else {
    fitMeasures$srmr <- NA
  }

  # information criteria:


  # Deviance based AIC / BIC / ebic (likelihood-based):
  fitMeasures <- add_ll_ic(fitMeasures)

  # Chi-square based AIC with df penalty (Kaplan, 2000): AIC(null) - AIC(saturated)
  fitMeasures$aic.x <- Tm - 2*fitMeasures$df

  # Chi-square based AIC with parameter penalty (Kline, 2016) - couldn't find the proper derivation
  fitMeasures$aic.x2 <- Tm + 2*fitMeasures$npar
  # fitMeasures$ebicTuning <- ebicTuning

  # Put in objet:
  x@fitmeasures <- fitMeasures
  return(x)
}
# 
# print.precisionFit <- function(x,...){
#   name <- deparse(substitute(x))[[1]]
#   if (nchar(name) > 10) name <- "object"
#   if (name=="x") name <- "object"
#   
#   cat("\nprecisionFit object:\n",
#       paste0("Use plot(",name,") to plot the network structure"),
#       "\n",
#       paste0("Fit measures stored under ",name,"$fitMeasures"),
#       "\n\n"
#   )
#   
#   fit <- data.frame(Measure = names(x$fitMeasures),
#                     Value = goodNum(unlist(x$fitMeasures)))
#   rownames(fit) <- NULL
#   print(fit)
# }
# 
# plot.precisionFit <- function(x,...){
#   qgraph::qgraph(x$network,...)
# }
