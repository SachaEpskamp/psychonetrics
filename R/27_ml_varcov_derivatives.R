# Analytic Jacobians of the distribution parameters phi with respect to the
# model parameters theta for ml_varcov models. Two variants exist, one per
# estimator (mirroring d_phi_theta_ml_lvm):
#
# - Two-level sufficient-statistics ML (estimator "ML"):
#     phi = [ mu ; vech(Sigma_within) ; vech(Sigma_between) ]
# - Wide-format FIML (estimator "FIML"):
#     phi = [ wide mu ; vech(wide Sigma) ], with
#     wide Sigma = I_nMax (x) Sigma_within + J_nMax (x) Sigma_between subset by
#     the design; the Jacobian is assembled in the compact row space
#     [ mu (p) ; vech(within-position block) (kp) ; vec(between-positions
#     block) (p^2) ] and mapped to the wide rows by the permutation matrix P
#     built in the constructor.
#
# In both variants the mean block is the identity and each covariance block is
# the single-level varcov covariance-structure Jacobian (the shared d_sigma_*
# helpers) for that level's type. There is no lambda, no beta, no residual and
# no latent-mean block -- ml_varcov carries only mu and the two covariance
# parameterizations. For every type the covariance block contributes p(p+1)/2
# parameters (for ggm that is p(p-1)/2 omega off-diagonals + p delta diagonals,
# and likewise for cor).

# Covariance-structure augmentation for one block, given its type and the
# helper matrices stored by impliedcovstructures under the block's suffix.
# Returns a (p(p+1)/2) x (p(p+1)/2) matrix mapping the block's covariance
# parameters to vech(Sigma_<block>). With lambda = I the d_sigma_sigma_zeta
# factor of the ml_lvm equivalent is the identity, so only the varcov
# augmentation survives:
aug_ml_varcov_block <- function(type, kp, L, In, C, D, A, Dstar,
                                sigma, lowertri, delta, SD,
                                delta_IminOinv, SD_IplusRho){
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

# Two-level (estimator "ML") group Jacobian:
#   rows    [ mu (p) ; vech Sigma_within (kp) ; vech Sigma_between (kp) ]
#   columns [ mu (p) ; within params (kp)     ; between params (kp)     ]
d_phi_theta_ml_varcov_group <- function(within, between, ...){
  dots <- list(...)

  # Number of variables (from the p x p within block; under FIML dots$mu is
  # the wide mean vector, so mu cannot be used for dimensions):
  nVar <- nrow(as.matrix(dots$sigma_within))
  kp <- nVar * (nVar + 1) / 2

  aug_within <- aug_ml_varcov_block(within, kp,
                                    L = dots$L, In = dots$In, C = dots$C, D = dots$D,
                                    A = dots$A, Dstar = dots$Dstar,
                                    sigma          = dots$sigma_within,
                                    lowertri       = dots$lowertri_within,
                                    delta          = dots$delta_within,
                                    SD             = dots$SD_within,
                                    delta_IminOinv = dots$delta_IminOinv_within,
                                    SD_IplusRho    = dots$SD_IplusRho_within)

  aug_between <- aug_ml_varcov_block(between, kp,
                                     L = dots$L, In = dots$In, C = dots$C, D = dots$D,
                                     A = dots$A, Dstar = dots$Dstar,
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

# Wide-format (estimator "FIML") group Jacobian. The compact rows are
#   [ mu (p) ; vech(within-position block) (kp) ; vec(between-positions block)
#     (p^2) ]
# where the within-position (diagonal) block of the wide covariance matrix is
# Sigma_within + Sigma_between and every between-positions (off-diagonal) block
# is Sigma_between. P (from the constructor) maps these compact rows to the
# wide distribution parameters [wide mu; vech(wide Sigma)]:
d_phi_theta_ml_varcov_wide_group <- function(within, between, ...){
  dots <- list(...)
  P <- dots$P

  nVar <- nrow(as.matrix(dots$sigma_within))
  kp <- nVar * (nVar + 1) / 2

  aug_within <- aug_ml_varcov_block(within, kp,
                                    L = dots$L, In = dots$In, C = dots$C, D = dots$D,
                                    A = dots$A, Dstar = dots$Dstar,
                                    sigma          = dots$sigma_within,
                                    lowertri       = dots$lowertri_within,
                                    delta          = dots$delta_within,
                                    SD             = dots$SD_within,
                                    delta_IminOinv = dots$delta_IminOinv_within,
                                    SD_IplusRho    = dots$SD_IplusRho_within)

  aug_between <- aug_ml_varcov_block(between, kp,
                                     L = dots$L, In = dots$In, C = dots$C, D = dots$D,
                                     A = dots$A, Dstar = dots$Dstar,
                                     sigma          = dots$sigma_between,
                                     lowertri       = dots$lowertri_between,
                                     delta          = dots$delta_between,
                                     SD             = dots$SD_between,
                                     delta_IminOinv = dots$delta_IminOinv_between,
                                     SD_IplusRho    = dots$SD_IplusRho_between)

  # Compact Jacobian:
  nobs     <- nVar + kp + nVar^2
  nelement <- nVar + kp + kp
  Jac <- Matrix(0, nrow = nobs, ncol = nelement, sparse = FALSE)

  meanInds       <- seq_len(nVar)
  inPositionInds <- nVar + seq_len(kp)
  betweenPosInds <- nVar + kp + seq_len(nVar^2)

  muInds      <- seq_len(nVar)
  withinInds  <- nVar + seq_len(kp)
  betweenInds <- nVar + kp + seq_len(kp)

  # Mean part:
  Jac[meanInds, muInds] <- Diagonal(nVar)

  # Within-position (diagonal) blocks equal Sigma_within + Sigma_between, so
  # both parameter blocks contribute:
  Jac[inPositionInds, withinInds]  <- aug_within
  Jac[inPositionInds, betweenInds] <- aug_between

  # Between-positions (off-diagonal) blocks equal Sigma_between; the compact
  # rows are the full vec, so map vech -> vec with the duplication matrix:
  Jac[betweenPosInds, betweenInds] <- dots$D %*% aug_between

  # Map the compact rows to the wide distribution parameters:
  Jac <- P %*% Jac

  sparseordense(Jac)
}

# Full Jacobian across groups (block-diagonal), mirroring d_phi_theta_ml_lvm.
# The prepare function flags which estimator variant applies; if the flag is
# absent (e.g. a C++-prepared prep), the wide variant is recognized by the P
# permutation matrix, which only FIML models carry:
d_phi_theta_ml_varcov <- function(prep){
  twolevel_ML <- if (!is.null(prep$twolevel_ML)){
    isTRUE(prep$twolevel_ML)
  } else {
    is.null(prep$groupModels[[1]]$P)
  }
  whatFun <- if (twolevel_ML){
    d_phi_theta_ml_varcov_group
  } else {
    d_phi_theta_ml_varcov_wide_group
  }

  d_per_group <- lapply(prep$groupModels, do.call, what = whatFun)

  if (length(d_per_group) == 1L){
    d_per_group[[1]]
  } else {
    bdiag(d_per_group)
  }
}
