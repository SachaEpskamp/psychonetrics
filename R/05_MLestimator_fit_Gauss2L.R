# Two-level (random intercept) sufficient-statistics ML estimator for
# multivariate normal data with cluster structure (ml_lvm, estimator = "ML").
#
# Implemented clean-room from the published two-level likelihood
# decomposition (McDonald & Goldstein, 1989; Muthen, 1990):
#
#   -2 l* = (N - J) [ ln|Sigma_W| + tr(Sigma_W^-1 S_PW) ] +
#           sum_s m_s [ ln|Sigma_s| + n_s tr(Sigma_s^-1 A_s) ],
#
# with Sigma_s = Sigma_W + n_s Sigma_B, A_s = cov.d_s +
# (mean.d_s - mu)(mean.d_s - mu)', S_PW the pooled within-cluster covariance
# (ML denominator N - J), and per distinct cluster size n_s: m_s clusters,
# mean.d_s the average of their cluster means and cov.d_s the ML covariance of
# those cluster means. The 2*pi constant is omitted here so that this fit
# function is numerically identical to the existing wide-format FIML objective
# (fimlEstimator_Gauss) at identical parameter values; see
# logLikelihood_gaussian2L_group for the full log-likelihood.

# -2 l* for one group, WITHOUT the N*p*log(2*pi) constant:
minustwo_logl_Gauss2L_noconstant <- function(mu, sigma_within, sigma_between, twolevel){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)

  N <- twolevel$N
  J <- twolevel$J

  # Non positive-definite within covariance: return a large penalty, mirroring
  # maxLikEstimator_Gauss_group:
  if (!sympd_cpp(SW)){
    return(1e20)
  }

  # Pooled-within part (zero weight when all clusters have size 1):
  res <- 0
  if (N > J){
    iSW <- solve_symmetric(SW)
    res <- (N - J) * (as.numeric(determinant(SW, logarithm = TRUE)$modulus) + sum(iSW * twolevel$S_PW))
  }

  # Between part, per distinct cluster size:
  sizes <- twolevel$sizes
  for (s in seq_len(nrow(sizes))){
    nj <- sizes$nj[s]
    m  <- sizes$m[s]
    Sj <- SW + nj * SB
    if (!sympd_cpp(Sj)){
      return(1e20)
    }
    iSj <- solve_symmetric(Sj)
    yc <- twolevel$mean_d[[s]] - mu
    A  <- twolevel$cov_d[[s]] + outer(yc, yc)
    res <- res + m * (as.numeric(determinant(Sj, logarithm = TRUE)$modulus) + nj * sum(iSj * A))
  }

  as.numeric(res)
}

# Fit function per group: (-2 l*) / J (no 2*pi constant):
maxLikEstimator_Gauss2L_group <- function(mu, sigma_within, sigma_between, twolevel, ...){
  minustwo_logl_Gauss2L_noconstant(mu = mu, sigma_within = sigma_within,
                                   sigma_between = sigma_between, twolevel = twolevel) / twolevel$J
}

# Fit function for the two-level Gaussian ML estimator. Total fit is the
# nobs-weighted (= clusters-weighted) average of the per-group fits, exactly
# as for the other estimators:
maxLikEstimator_Gauss2L <- function(prep){
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels, do.call, what = maxLikEstimator_Gauss2L_group)
  sum(fit_per_group)
}
