# Implied model for ml_varcov. ml_varcov is the multi-level variance-covariance
# model: the observed within- and between-cluster covariance matrices are
# modelled directly (no latent layer). With the observed decomposition
# y_ij = mu + eta^W_ij + eta^B_i,
#   var(y_ij)            = sigma_within + sigma_between
#   cov(y_ij, y_ij')     = sigma_between      (j != j', same cluster)
# and each of sigma_within / sigma_between is parameterized by its own type
# (cov / chol / prec / ggm / cor) exactly as in the single-level varcov family.
#
# What is produced depends on the estimator:
# - estimator "ML" (distribution "TwoLevelGaussian"): the two-level
#   distribution parameters mu (p), sigma_within (p x p), sigma_between (p x p)
#   consumed by the Gauss2L estimator.
# - estimator "FIML" (distribution "Gaussian"): additionally the wide-format
#   (one row per cluster) mean vector and covariance matrix
#     fullSigma = I_nMax (x) sigma_within + J_nMax (x) sigma_between
#   subset by the design pattern, plus its inverse, stored under mu / sigma /
#   kappa for the per-pattern FIML estimator (mirrors the ml_lvm wide format).
implied_ml_varcov <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)
  }

  # Implied covariance structures for the two blocks. impliedcovstructures reads
  # the {omega,delta,kappa,lowertri,rho,SD}_<name> matrices and writes
  # sigma_<name> (plus the derivative helpers delta_IminOinv_<name> etc. when
  # all = FALSE), so after these two calls x[[g]] already holds sigma_within and
  # sigma_between under their native names:
  if (model@cpp){
    x <- impliedcovstructures_cpp(x, "within",  type = model@types$within,  all = all)
    x <- impliedcovstructures_cpp(x, "between", type = model@types$between, all = all)
  } else {
    x <- impliedcovstructures(x, "within",  type = model@types$within,  all = all)
    x <- impliedcovstructures(x, "between", type = model@types$between, all = all)
  }

  # Wide format needed for the FIML estimator:
  wide <- model@estimator == "FIML"
  if (wide){
    designPattern <- model@extramatrices$designPattern
    nMaxInCluster <- ncol(designPattern)
  }

  for (g in seq_along(x)){
    # Two-level distribution parameters (force symmetric, mirroring ml_lvm/ml_var1):
    mu_p <- as.matrix(x[[g]]$mu)
    sw <- as.matrix(x[[g]]$sigma_within)
    sb <- as.matrix(x[[g]]$sigma_between)
    sw <- 0.5 * (sw + t(sw))
    sb <- 0.5 * (sb + t(sb))
    x[[g]]$mu <- mu_p
    x[[g]]$sigma_within  <- sw
    x[[g]]$sigma_between <- sb

    if (wide){
      # Wide-format mean vector and covariance matrix, subset by the design
      # (also under all = TRUE, so that getmatrix() returns the wide sigma /
      # kappa / mu exactly as it does for ml_lvm FIML models):
      sub <- as.vector(designPattern) == 1
      fullMu <- do.call(rbind, lapply(seq_len(nMaxInCluster), function(t){
        mu_p[designPattern[, t] == 1, , drop = FALSE]
      }))

      fullSigma <- as.matrix(Diagonal(nMaxInCluster) %x% sw +
                               matrix(1, nMaxInCluster, nMaxInCluster) %x% sb)
      fullSigma <- fullSigma[sub, sub, drop = FALSE]
      fullSigma <- 0.5 * (fullSigma + t(fullSigma))

      x[[g]]$mu <- as.matrix(fullMu)
      x[[g]]$sigma <- fullSigma
      x[[g]]$kappa <- solve_symmetric(fullSigma, logdet = TRUE)
    }

    if (all){
      # Cross-sectional (total) covariance is the sum of the two levels:
      x[[g]]$sigma_crosssection <- x[[g]]$sigma_within + x[[g]]$sigma_between
    }
  }

  x
}
