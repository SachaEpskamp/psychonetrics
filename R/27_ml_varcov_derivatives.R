# Analytic Jacobian of the two-level distribution parameters
#   phi = [ mu ; vech(Sigma_within) ; vech(Sigma_between) ]
# with respect to the model parameters theta, for a single group of an
# ml_varcov model. This is the multi-level analogue of d_phi_theta_varcov_group:
# the mean block is the identity, and each covariance block is the single-level
# varcov covariance-structure Jacobian (the shared d_sigma_* helpers) for that
# level's type. There is no lambda, no beta, no residual and no latent-mean
# block -- ml_varcov carries only mu and the two covariance parameterizations.
#
# Row order (nobs = p + p(p+1)/2 + p(p+1)/2):
#   [ mu ; vech(Sigma_within) ; vech(Sigma_between) ]
# Column (theta) order, fixed by the matrix registration order in the
# constructor (mu, then the "within" block, then the "between" block):
#   [ mu ; within covariance parameters ; between covariance parameters ]
# For every type the covariance block contributes p(p+1)/2 parameters (for ggm
# that is p(p-1)/2 omega off-diagonals + p delta diagonals, and likewise for
# cor), so the within and between column blocks are each p(p+1)/2 wide.
d_phi_theta_ml_varcov_group <- function(within, between, ...){
  dots <- list(...)

  # Number of variables (mu is p-variate):
  nVar <- length(as.vector(dots$mu))
  kp <- nVar * (nVar + 1) / 2

  # Shared p-dimensional derivative kernels (from model@extramatrices):
  L     <- dots$L      # elimination matrix (vec -> vech)
  In    <- dots$In     # identity p x p
  A     <- dots$A      # diagonalization matrix
  Dstar <- dots$Dstar  # strict duplication matrix (off-diagonal only)
  D     <- dots$D      # duplication matrix (vech -> vec)
  C     <- dots$C      # commutation matrix

  # Covariance-structure augmentation for one block, given its type and the
  # helper matrices stored by impliedcovstructures under the block's suffix.
  # Returns a (p(p+1)/2) x (p(p+1)/2) matrix mapping the block's covariance
  # parameters to vech(Sigma_<block>). Mirrors the aug blocks of
  # d_phi_theta_ml_lvm2L_group (with lambda = I the d_sigma_sigma_zeta factor is
  # the identity, so only the varcov augmentation survives).
  aug <- function(type, sigma, lowertri, delta, SD, delta_IminOinv, SD_IplusRho){
    if (type == "cov"){
      Diagonal(kp)
    } else if (type == "chol"){
      d_sigma_cholesky(lowertri = lowertri, L = L, C = C, In = In)
    } else if (type == "prec"){
      d_sigma_kappa(L = L, D = D, sigma = sigma)
    } else if (type == "ggm"){
      cbind(
        d_sigma_omega(L = L, delta_IminOinv = delta_IminOinv, A = A, delta = delta, Dstar = Dstar),
        d_sigma_delta(L = L, delta_IminOinv = delta_IminOinv, In = In, delta = delta, A = A)
      )
    } else if (type == "cor"){
      cbind(
        d_sigma_rho(L = L, SD = SD, A = A, Dstar = Dstar),
        d_sigma_SD(L = L, SD_IplusRho = SD_IplusRho, In = In, A = A)
      )
    } else {
      stop("Unsupported ml_varcov type: ", type)
    }
  }

  aug_within <- aug(within,
                    sigma          = dots$sigma_within,
                    lowertri       = dots$lowertri_within,
                    delta          = dots$delta_within,
                    SD             = dots$SD_within,
                    delta_IminOinv = dots$delta_IminOinv_within,
                    SD_IplusRho    = dots$SD_IplusRho_within)

  aug_between <- aug(between,
                     sigma          = dots$sigma_between,
                     lowertri       = dots$lowertri_between,
                     delta          = dots$delta_between,
                     SD             = dots$SD_between,
                     delta_IminOinv = dots$delta_IminOinv_between,
                     SD_IplusRho    = dots$SD_IplusRho_between)

  # Assemble the block Jacobian:
  nobs     <- nVar + kp + kp
  nelement <- nVar + kp + kp
  Jac <- Matrix(0, nrow = nobs, ncol = nelement, sparse = FALSE)

  meanInds <- seq_len(nVar)
  swInds   <- nVar + seq_len(kp)
  sbInds   <- nVar + kp + seq_len(kp)

  muInds      <- seq_len(nVar)
  withinInds  <- nVar + seq_len(kp)
  betweenInds <- nVar + kp + seq_len(kp)

  # Mean part: mu enters phi's mean block directly (identity), and enters
  # neither covariance block:
  Jac[meanInds, muInds] <- Diagonal(nVar)

  # Within covariance block:
  Jac[swInds, withinInds] <- aug_within

  # Between covariance block:
  Jac[sbInds, betweenInds] <- aug_between

  sparseordense(Jac)
}

# Full Jacobian across groups (block-diagonal), mirroring d_phi_theta_ml_lvm:
d_phi_theta_ml_varcov <- function(prep){
  d_per_group <- lapply(prep$groupModels, do.call, what = d_phi_theta_ml_varcov_group)

  if (length(d_per_group) == 1L){
    d_per_group[[1]]
  } else {
    bdiag(d_per_group)
  }
}
