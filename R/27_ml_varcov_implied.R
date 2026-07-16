# Implied model for ml_varcov. Produces, per group, the two-level distribution
# parameters consumed by the Gauss2L estimator: mu (p), sigma_within (p x p),
# sigma_between (p x p). R-only (model@cpp is forced FALSE by the constructor
# in v1; the C++ twin implied_ml_varcov_cpp is registered for model@cpp = TRUE).
#
# ml_varcov is the multi-level variance-covariance model: the observed within-
# and between-cluster covariance matrices are modelled directly (no latent
# layer). With the observed decomposition y_ij = mu + eta^W_ij + eta^B_i,
#   var(y_ij)            = sigma_within + sigma_between
#   cov(y_ij, y_ij')     = sigma_between      (j != j', same cluster)
# and each of sigma_within / sigma_between is parameterized by its own type
# (cov / chol / prec / ggm / cor) exactly as in the single-level varcov family.
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

  for (g in seq_along(x)){
    # Two-level distribution parameters (force symmetric, mirroring ml_lvm/ml_var1):
    x[[g]]$mu <- as.matrix(x[[g]]$mu)
    sw <- as.matrix(x[[g]]$sigma_within)
    sb <- as.matrix(x[[g]]$sigma_between)
    x[[g]]$sigma_within  <- 0.5 * (sw + t(sw))
    x[[g]]$sigma_between <- 0.5 * (sb + t(sb))

    if (all){
      # Cross-sectional (total) covariance is the sum of the two levels:
      x[[g]]$sigma_crosssection <- x[[g]]$sigma_within + x[[g]]$sigma_between
    }
  }

  x
}
