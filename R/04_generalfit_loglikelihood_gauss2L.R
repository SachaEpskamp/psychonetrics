# Log-likelihood for the two-level Gaussian distribution
# (ml_lvm with estimator = "ML"), including the 2*pi constant:
#
#   loglik = -0.5 * ( -2 l*  +  (#observed values) * log(2*pi) ),
#
# with -2 l* the kernel computed in minustwo_logl_Gauss2L_noconstant
# (05_MLestimator_fit_Gauss2L.R). For COMPLETE data the number of observed
# values is N * p (N level-1 units, p variables); for MISSING data it is the
# total number of non-NA scalar values, stored as twolevel$nel.

# Full log-likelihood per group:
logLikelihood_gaussian2L_group <- function(mu, sigma_within, sigma_between, twolevel, ...){
  p <- length(as.vector(mu))
  m2ll <- minustwo_logl_Gauss2L_noconstant(mu = mu, sigma_within = sigma_within,
                                           sigma_between = sigma_between, twolevel = twolevel)
  nel <- if (isTRUE(twolevel$missing)) twolevel$nel else twolevel$N * p
  -0.5 * (m2ll + nel * log(2*pi))
}

# Log-likelihood over all groups:
logLikelihood_gaussian2L <- function(prep){
  ll_per_group <- sapply(prep$groupModels, do.call, what = logLikelihood_gaussian2L_group)
  sum(ll_per_group)
}
