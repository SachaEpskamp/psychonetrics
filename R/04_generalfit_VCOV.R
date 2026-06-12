getVCOV <- function(model,approximate_SEs = FALSE){
  # The two-level ML estimator only has an R implementation:
  model <- force_R_path_if_needed(model)

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
