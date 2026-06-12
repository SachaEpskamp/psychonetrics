# Log-likelihood for the two-level Gaussian distribution
# (ml_lvm with estimator = "ML"), including the 2*pi constant:
#
#   loglik = -0.5 * ( -2 l*  +  N * p * log(2*pi) ),
#
# with -2 l* the kernel computed in minustwo_logl_Gauss2L_noconstant
# (05_MLestimator_fit_Gauss2L.R) and N the total number of level-1 units.

# Full log-likelihood per group:
logLikelihood_gaussian2L_group <- function(mu, sigma_within, sigma_between, twolevel, ...){
  p <- length(as.vector(mu))
  m2ll <- minustwo_logl_Gauss2L_noconstant(mu = mu, sigma_within = sigma_within,
                                           sigma_between = sigma_between, twolevel = twolevel)
  -0.5 * (m2ll + twolevel$N * p * log(2*pi))
}

# Log-likelihood over all groups:
logLikelihood_gaussian2L <- function(prep){
  ll_per_group <- sapply(prep$groupModels, do.call, what = logLikelihood_gaussian2L_group)
  sum(ll_per_group)
}
