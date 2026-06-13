# Analytic expected Hessian of the two-level Gaussian ML fit function (see
# 05_MLestimator_fit_Gauss2L.R) with respect to the distribution parameters
# phi = [mu; vech(Sigma_W); vech(Sigma_B)]. Clean-room implementation from the
# expected per-cluster-unit Fisher information of the two-level likelihood
# (verified against lavaan's lav_mvnorm_cluster_information_expected kernel).
#
# With Omega_s = Sigma_B + Sigma_W / n_s = (Sigma_W + n_s Sigma_B) / n_s the
# expected information per cluster is
#
#   I(mu)          = (1/J) sum_s m_s Omega_s^-1
#   I(vech blocks) = (1/J) [ sum_s m_s D_s' (0.5 D'(Omega_s^-1 (x) Omega_s^-1) D) D_s
#                    + (N - J) blkdiag( 0.5 D'(Sigma_W^-1 (x) Sigma_W^-1) D, 0 ) ]
#
# with D the duplication matrix, k = p(p+1)/2 and D_s = [ (1/n_s) I_k , I_k ]
# (the size-dependent map vech Omega_s = (1/n_s) vech Sigma_W + vech Sigma_B).
# The mean block is orthogonal to the covariance blocks.
#
# Following the convention of the other Gaussian expected-Hessian functions
# (e.g. expected_hessian_Gaussian_group: 2*kappa for the means and
# D'(kappa (x) kappa)D for the (co)variances, both equal to TWICE the per-unit
# expected information), this function returns 2 * I per group, so that the
# generic assembly 0.5 * M'J'HJM in psychonetrics_FisherInformation yields the
# per-cluster unit information and getVCOV's 1/n scaling (n = total number of
# clusters) gives the correct asymptotic covariance matrix.

# Per group:
expected_hessian_Gauss2L_group <- function(mu, sigma_within, sigma_between, twolevel, D_y, ...){
  # The closed-form expected information below assumes COMPLETE data within
  # each distinct cluster size. With within-cluster missingness the expected
  # information is pattern-dependent and has no such closed form; those models
  # therefore use numeric Fisher information for SEs (see
  # psychonetrics_FisherInformation / runmodel). This branch should not be
  # reached for missing data:
  if (isTRUE(twolevel$missing)){
    stop("Internal error: analytic expected Hessian is not available for two-level ML with missing data; numeric Fisher information should be used instead. Please report this bug.")
  }
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)
  p <- length(mu)
  k <- p * (p + 1) / 2

  N <- twolevel$N
  J <- twolevel$J

  # Duplication matrix for the p observed variables:
  D <- D_y

  # Mean part and the three covariance blocks (WW, WB, BB):
  I_mu <- matrix(0, p, p)
  I_WW <- matrix(0, k, k)
  I_WB <- matrix(0, k, k)
  I_BB <- matrix(0, k, k)

  sizes <- twolevel$sizes
  for (s in seq_len(nrow(sizes))){
    nj <- sizes$nj[s]
    m  <- sizes$m[s]

    # Omega_s^-1 = n_s * (Sigma_W + n_s Sigma_B)^-1:
    iSj <- solve_symmetric(SW + nj * SB)
    iOm <- nj * iSj

    # Mean part: m_s * Omega_s^-1:
    I_mu <- I_mu + m * iOm

    # G_s = 0.5 D'(Sigma_s^-1 (x) Sigma_s^-1)D; with D_s = [(1/n_s) I, I] and
    # H_s = 0.5 D'(Omega_s^-1 (x) Omega_s^-1)D = n_s^2 G_s, the blocks of
    # D_s' H_s D_s are: WW = G_s, WB = n_s G_s, BB = n_s^2 G_s:
    G <- as.matrix(0.5 * t(D) %*% (iSj %x% iSj) %*% D)

    I_WW <- I_WW + m * G
    I_WB <- I_WB + m * nj * G
    I_BB <- I_BB + m * nj^2 * G
  }

  # Pooled-within part (zero weight when all clusters have size 1):
  if (N > J){
    iSW <- solve_symmetric(SW)
    I_WW <- I_WW + (N - J) * as.matrix(0.5 * t(D) %*% (iSW %x% iSW) %*% D)
  }

  # Assemble the covariance part:
  I_cov <- rbind(
    cbind(I_WW, I_WB),
    cbind(I_WB, I_BB)
  )

  # Return TWICE the per-cluster unit information (see header):
  bdiag(2 / J * I_mu, 2 / J * I_cov)
}

# Total (mirrors expected_hessian_Gaussian):
expected_hessian_Gauss2L <- function(prep){
  # model is already prepared!

  # expected Hessian per group:
  exph_per_group <- lapply(prep$groupModels, do.call, what = expected_hessian_Gauss2L_group)

  # Weight:
  for (i in seq_along(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }

  # Bind by block and return:
  Reduce("bdiag", exph_per_group)
}
