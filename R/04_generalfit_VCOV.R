# Robust ML standard errors (Phase 1, complete data). Returns the parameter
# covariance matrix (already divided by N) for a robust ML estimator, or NULL
# if the required building blocks are unavailable (so the caller falls back to
# the naive VCOV). The point estimates are plain ML; only the covariance of the
# estimator changes.
#
#  robust.sem (Browne 1984; MLM/MLMV/MLMVS):
#     Vcov = E^-1 [ sum_g f_g Delta_g' V_g Gamma_g V_g Delta_g ] E^-1 / N
#  robust.huber.white (MLR):
#     Vcov = A^-1 [ sum_g f_g Delta_g' B1_g Delta_g ] A^-1 / N
#  with A the OBSERVED unit information (or x@information when
#  information = "expected").
getVCOV_robust <- function(model, information = c("observed","expected")){
  information <- match.arg(information)
  cfg <- get_robust_config(model)
  se_type <- cfg$se
  if (is.null(se_type) || !nzchar(se_type)) return(NULL)

  n <- sum(model@sample@groups$nobs)

  if (se_type == "robust.sem"){
    comp <- tryCatch(ml_robust_components(model), error = function(e) NULL)
    if (is.null(comp)) return(NULL)
    E_inv <- tryCatch(solve(comp$E), error = function(e) NULL)
    if (is.null(E_inv) || any(!is.finite(E_inv))) return(NULL)
    meat <- matrix(0, comp$npar, comp$npar)
    for (g in seq_len(comp$nGroups)){
      VD <- comp$V[[g]] %*% comp$Delta[[g]]
      meat <- meat + comp$fg[g] * (t(VD) %*% comp$Gamma[[g]] %*% VD)
    }
    return(E_inv %*% meat %*% E_inv / n)
  }

  if (se_type == "robust.huber.white"){
    meat <- tryCatch(mlr_meat_components(model), error = function(e) NULL)
    if (is.null(meat)) return(NULL)
    # Bread A: observed information (hessian-based) by default; the expected
    # information (x@information) is offered as a cheaper alternative.
    if (information == "expected"){
      if (!is.null(model@information) && length(model@information) > 0){
        A <- as.matrix(model@information)
      } else if (model@cpp){
        A <- psychonetrics_FisherInformation_cpp(model)
      } else {
        A <- psychonetrics_FisherInformation(model)
      }
    } else {
      A <- meat$A
    }
    A_inv <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(A_inv) || any(!is.finite(A_inv))) return(NULL)
    B0 <- matrix(0, meat$npar, meat$npar)
    for (g in seq_len(meat$nGroups)){
      B0 <- B0 + meat$fg[g] * (t(meat$Delta[[g]]) %*% meat$B1[[g]] %*% meat$Delta[[g]])
    }
    return(A_inv %*% B0 %*% A_inv / n)
  }

  NULL
}

getVCOV <- function(model,approximate_SEs = FALSE){

  # Robust ML standard errors (MLM/MLMV/MLMVS/MLR). These map internally to
  # estimator = "ML" (complete data) or, for MLR with within-row missing data,
  # estimator = "FIML"; the point estimates are unchanged but the sampling
  # covariance is a sandwich estimator. Fall through to the naive VCOV below if
  # the robust building blocks are unavailable.
  if (model@estimator %in% c("ML", "FIML") && is_robust_ML(model)){
    robInfo <- tryCatch(getVCOV_robust(model), error = function(e) NULL)
    if (!is.null(robInfo) && all(is.finite(robInfo))){
      return(robInfo)
    }
    # else: fall back to naive ML VCOV (computed below).
  }

  # Wishart Gaussian likelihood (varcov / lvm, complete-data ML): the parameter
  # covariance is scaled by 1/(n_g - 1) per group instead of 1/n_g (so SEs are
  # inflated by sqrt(n_g/(n_g-1))), matching lavaan likelihood = "wishart". The
  # unit Fisher information weights group g by n_g/N (see expected_hessian_*),
  # so recomputing it on a copy of the model whose group sizes are reduced to
  # (n_g - 1) yields sum_g ((n_g-1)/(N-G)) I_g; dividing that by (N-G) gives
  # exactly (sum_g (n_g-1) I_g)^-1 — the Wishart sampling covariance for any
  # number of groups. The point estimates are unaffected (the implied matrices
  # do not depend on nobs).
  if (is_wishart(model) && model@estimator %in% c("ML")){
    nobs_g <- model@sample@groups$nobs
    if (all(nobs_g > 1)){
      model_w <- model
      model_w@sample@groups$nobs <- nobs_g - 1
      # Force a fresh information computation on the reduced-nobs copy (the
      # cached x@information used the full nobs):
      model_w@information <- matrix(0, 0, 0)
      if (model_w@cpp){
        Info_w <- psychonetrics_FisherInformation_cpp(model_w)
      } else {
        Info_w <- psychonetrics_FisherInformation(model_w)
      }
      nW <- sum(nobs_g - 1)
      return(1/nW * as.matrix(solve_symmetric(Info_w, approx = approximate_SEs)))
    }
  }

  # The @information slot is a "matrix" slot, so when a model is run with
  # addInformation = FALSE it holds the default empty 0x0 matrix rather than
  # NULL. Treat a zero-size matrix as absent and recompute the Fisher
  # information on demand (mirrors addSEs_cpp()).
  if (!is.null(model@information) && length(model@information) > 0){
    Info <- model@information
  } else {
    if (model@cpp){
      Info <- psychonetrics_FisherInformation_cpp(model)
    } else {
      Info <- psychonetrics_FisherInformation(model)
    }

  }

  # Total sample size:
  n <- sum(model@sample@groups$nobs)

  # Naive (model-based) VCOV: (1/n) * Info^-1. For ML/FIML/WLS this is the
  # correct asymptotic covariance. For full WLS the weight matrix equals
  # Gamma^-1, so Info = Delta'Gamma^-1 Delta and the naive form is already the
  # sandwich.
  naive <- 1/n * as.matrix(solve_symmetric(Info,approx=approximate_SEs))

  # For the least-squares estimators ULS and DWLS the weight matrix W used in
  # estimation is NOT equal to Gamma^-1, so the naive (Delta'WDelta)^-1 form is
  # not a consistent estimator of the sampling covariance. Replace it with the
  # robust sandwich
  #   (1/n) (Delta'WDelta)^-1 (Delta'W Gamma W Delta) (Delta'WDelta)^-1
  # which reduces to the naive form when W = Gamma^-1 (full WLS). See
  # wls_sandwich_components() for the per-group Delta/W/Gamma assembly (shared
  # with the WLSMV scaled test statistic), which already applies the n_g/n
  # group weighting consistent with Info.
  if (model@estimator %in% c("ULS","DWLS")){
    comp <- tryCatch(wls_sandwich_components(model), error = function(e) NULL)

    if (is.null(comp)){
      # Gamma unavailable (e.g. a model fit to a covariance matrix rather than
      # raw data): cannot form the sandwich. Fall back to the naive VCOV but
      # warn that the standard errors are not robust.
      warning("Robust (sandwich) standard errors for the ", model@estimator,
              " estimator require the asymptotic covariance matrix (Gamma) of ",
              "the sample statistics, which is unavailable (e.g. when a model ",
              "is fit to a covariance matrix rather than to raw data). ",
              "Returning non-robust standard errors, which assume the weight ",
              "matrix equals the inverse of Gamma and are generally inaccurate ",
              "for ULS/DWLS.")
      return(naive)
    }

    # Bread inverse E_inv = (sum_g Delta_g' W_g Delta_g)^-1:
    E_inv <- tryCatch(solve(comp$DtWD), error = function(e) NULL)
    if (is.null(E_inv) || any(!is.finite(E_inv))){
      return(naive)
    }

    # Meat B = sum_g Delta_g' W_g Gamma_g W_g Delta_g (group-weighted):
    B <- matrix(0, comp$npar, comp$npar)
    for (g in seq_len(comp$nGroups)){
      WD <- comp$W[[g]] %*% comp$Delta[[g]]
      B <- B + t(WD) %*% comp$Gamma[[g]] %*% WD
    }

    # Sandwich VCOV (note: comp$DtWD and B already carry the 1/n-consistent
    # group weighting, so the overall scale matches the naive (1/n) Info^-1):
    return((1/n) * E_inv %*% B %*% E_inv)
  }

  # Attempt to invert (ML / FIML / WLS):
  naive
}
