# Gradient of the two-level Gaussian ML fit function (see
# 05_MLestimator_fit_Gauss2L.R) with respect to the distribution parameters
# phi = [mu; vech(Sigma_W); vech(Sigma_B)]. Clean-room implementation from the
# same published formulas. With K_s = Sigma_s^-1 - n_s Sigma_s^-1 A_s Sigma_s^-1:
#
#   d(-2 l*)/d mu       = -2 sum_s m_s n_s Sigma_s^-1 (mean.d_s - mu)
#   d(-2 l*)/d Sigma_W  = sum_s m_s K_s + (N - J)(Sigma_W^-1 - Sigma_W^-1 S_PW Sigma_W^-1)
#   d(-2 l*)/d Sigma_B  = sum_s m_s n_s K_s
#
# Matrix derivatives map to vech by doubling the off-diagonal elements and
# keeping the diagonal (standard duplication-matrix convention, matching the
# vech/elimination ordering used elsewhere in psychonetrics).

# Helper: matrix derivative -> vech gradient (double off-diagonal entries):
vech_grad_2L <- function(G){
  G2 <- 2 * G
  diag(G2) <- diag(G)
  G2[lower.tri(G2, diag = TRUE)]
}

# Jacobian (1 x npar row vector) per group: (1/J) d(-2 l*)/d phi.
# The group weighting (J_g / J_total) is applied by the outer function, so
# together the gradient matches the fit function
# F_g = (-2 l*_g) / J_g, F = sum_g (J_g/J_total) F_g:
jacobian_gaussian2L_group_sigma <- function(mu, sigma_within, sigma_between, twolevel, ...){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)
  p <- length(mu)

  N <- twolevel$N
  J <- twolevel$J

  g_mu <- numeric(p)
  G_SW <- matrix(0, p, p)
  G_SB <- matrix(0, p, p)

  sizes <- twolevel$sizes
  for (s in seq_len(nrow(sizes))){
    nj <- sizes$nj[s]
    m  <- sizes$m[s]
    iSj <- solve_symmetric(SW + nj * SB)
    yc <- twolevel$mean_d[[s]] - mu
    A  <- twolevel$cov_d[[s]] + outer(yc, yc)
    K  <- iSj - nj * iSj %*% A %*% iSj
    g_mu <- g_mu - m * 2 * nj * as.numeric(iSj %*% yc)
    G_SW <- G_SW + m * K
    G_SB <- G_SB + m * nj * K
  }

  if (N > J){
    iSW <- solve_symmetric(SW)
    G_SW <- G_SW + (N - J) * (iSW - iSW %*% twolevel$S_PW %*% iSW)
  }

  grad <- c(g_mu, vech_grad_2L(G_SW), vech_grad_2L(G_SB)) / J

  # Return as 1 x npar row matrix (mirrors jacobian_gaussian_group_sigma):
  matrix(grad, nrow = 1)
}

# Now for all groups (mirrors jacobian_gaussian_sigma):
jacobian_gaussian2L_sigma <- function(prep){
  # Jacobian per group:
  g_per_group <- lapply(prep$groupModels, do.call, what = jacobian_gaussian2L_group_sigma)

  # Weight:
  for (i in seq_along(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }

  # Bind by column and return:
  Reduce("cbind", g_per_group)
}
