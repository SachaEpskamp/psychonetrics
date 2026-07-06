# Model Jacobian d phi / d theta for ml_var1, where per group
#   phi   = [ mu_z (2p) ; vech(Sigma_W) (k) ; vech(Sigma_Bf) (k) ],  k = p(2p+1)
#   theta = [ mu (p) ; vec(beta) (p^2) ; zeta_within block (kp) ; zeta_between (kp) ]
# with kp = p(p+1)/2 and vech = lower-triangle (incl. diagonal), column-major.
#
# The within block reuses the exported var1 derivative kernels
# (d_sigma0_beta_var1_cpp etc.), mapped from [vech(Sigma0); vec(Sigma1)] to
# vech(Sigma_W) via the P_within structure map; the between block maps
# vech(Sigma_B) to vech(1_{2x2} (x) Sigma_B) via P_between.
#
# In the free-within saturated mode (toeplitz = FALSE) the within/between
# blocks are free 2p Cholesky parameterizations and mu is a free 2p mean.

# Augmentation of a typed covariance block (returns d vech(sigma)/d params,
# with kp columns). 'suffix' selects the _zeta_within / _zeta_between helpers.
.ml_var1_augment <- function(type, dots, kp,
                             In, L, C, D, Dstar, A){
  if (type == "cov"){
    Diagonal(kp)
  } else if (type == "chol"){
    d_sigma_cholesky(lowertri = dots$lowertri, L = L, C = C, In = In)
  } else if (type == "prec"){
    d_sigma_kappa(L = L, D = D, sigma = dots$sigma)
  } else if (type == "ggm"){
    cbind(
      d_sigma_omega(L = L, delta_IminOinv = dots$delta_IminOinv, A = A, delta = dots$delta, Dstar = Dstar),
      d_sigma_delta(L = L, delta_IminOinv = dots$delta_IminOinv, In = In, delta = dots$delta, A = A)
    )
  } else if (type == "cor"){
    cbind(
      d_sigma_rho(L = L, SD = dots$SD, A = A, Dstar = Dstar),
      d_sigma_SD(L = L, SD_IplusRho = dots$SD_IplusRho, In = In, A = A)
    )
  }
}

d_phi_theta_ml_var1_group <- function(within_latent, between_latent, toeplitz, ...){
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Free-within saturated mode:
  # -------------------------------------------------------------------------
  if (!isTRUE(toeplitz)){
    L2p <- dots$L2p; C2p <- dots$C2p; I2p <- dots$I2p
    nVar <- nrow(I2p)                     # = 2p
    k <- nVar * (nVar + 1) / 2

    nobs <- nVar + k + k
    nelement <- nVar + k + k
    Jac <- Matrix(0, nrow = nobs, ncol = nelement, sparse = FALSE)

    meanInds <- seq_len(nVar)
    swInds <- nVar + seq_len(k)
    sbInds <- nVar + k + seq_len(k)

    muInds <- seq_len(nVar)
    cholWInds <- nVar + seq_len(k)
    cholBInds <- nVar + k + seq_len(k)

    Jac[meanInds, muInds] <- Diagonal(nVar)
    Jac[swInds, cholWInds] <- d_sigma_cholesky(lowertri = dots$lowertri_zeta_within, L = L2p, C = C2p, In = I2p)
    Jac[sbInds, cholBInds] <- d_sigma_cholesky(lowertri = dots$lowertri_zeta_between, L = L2p, C = C2p, In = I2p)

    return(sparseordense(Jac))
  }

  # -------------------------------------------------------------------------
  # Toeplitz (structured) mode:
  # -------------------------------------------------------------------------
  beta <- as.matrix(dots$beta)
  p <- nrow(beta)
  kp <- p * (p + 1) / 2
  k <- p * (2 * p + 1)

  In <- dots$In; L <- dots$L; C <- dots$C; D2 <- dots$D2; Dstar <- dots$Dstar; A <- dots$A
  P_within <- dots$P_within; P_between <- dots$P_between
  BetaStar <- dots$BetaStar; IkronBeta <- dots$IkronBeta
  sigma_within <- as.matrix(dots$sigma_within)

  nobs <- 2 * p + k + k
  nelement <- p + p^2 + kp + kp
  Jac <- Matrix(0, nrow = nobs, ncol = nelement, sparse = FALSE)

  # Row blocks (phi):
  meanInds <- seq_len(2 * p)
  swInds <- 2 * p + seq_len(k)
  sbInds <- 2 * p + k + seq_len(k)

  # Column blocks (theta):
  muInds  <- seq_len(p)
  betaInds <- p + seq_len(p^2)
  zwInds  <- p + p^2 + seq_len(kp)
  zbInds  <- p + p^2 + kp + seq_len(kp)

  # Augmentation matrices for the two typed blocks:
  aug_within <- .ml_var1_augment(within_latent, list(
    lowertri = dots$lowertri_zeta_within,
    sigma = dots$sigma_zeta_within,
    delta_IminOinv = dots$delta_IminOinv_zeta_within,
    delta = dots$delta_zeta_within,
    SD = dots$SD_zeta_within,
    SD_IplusRho = dots$SD_IplusRho_zeta_within
  ), kp, In, L, C, D2, Dstar, A)

  aug_between <- .ml_var1_augment(between_latent, list(
    lowertri = dots$lowertri_zeta_between,
    sigma = dots$sigma_zeta_between,
    delta_IminOinv = dots$delta_IminOinv_zeta_between,
    delta = dots$delta_zeta_between,
    SD = dots$SD_zeta_between,
    SD_IplusRho = dots$SD_IplusRho_zeta_between
  ), kp, In, L, C, D2, Dstar, A)

  ### Means: d mu_z / d mu = [I_p ; I_p] ###
  Jac[meanInds, muInds] <- rbind(Diagonal(p), Diagonal(p))

  ### Within block ###
  # d vech(Sigma0) / d vec(beta):
  Jb  <- d_sigma0_beta_var1_cpp(BetaStar = BetaStar, In = In, sigma = sigma_within, C = C, L = L)
  # d vech(Sigma0) / d (zeta_within params) (base matrix for the cpp kernels):
  J0z <- as.matrix(d_sigma0_sigma_zeta_var1_cpp(L = L, BetaStar = BetaStar, D2 = D2) %*% aug_within)
  # d vec(Sigma1) / d vec(beta):
  J1b <- d_sigma1_beta_var1_cpp(beta = beta, Jb = Jb, IkronBeta = IkronBeta, D2 = D2, sigma = sigma_within, In = In)
  # d vec(Sigma1) / d (zeta_within params) (pass the AUGMENTED J0z as Js):
  J1z <- d_sigma1_sigma_zeta_var1_cpp(Js = J0z, IkronBeta = IkronBeta, D2 = D2)

  # Map [vech(Sigma0); vec(Sigma1)] to vech(Sigma_W):
  Jac[swInds, betaInds] <- as.matrix(P_within %*% rbind(Jb, J1b))
  Jac[swInds, zwInds]   <- as.matrix(P_within %*% rbind(J0z, J1z))

  ### Between block: d vech(1_{2x2} (x) Sigma_B) / d (zeta_between params) ###
  Jac[sbInds, zbInds] <- as.matrix(P_between %*% aug_between)

  return(sparseordense(Jac))
}

# Full model Jacobian (block-diagonal across groups):
d_phi_theta_ml_var1 <- function(prep){
  d_per_group <- lapply(prep$groupModels, do.call, what = d_phi_theta_ml_var1_group)
  if (length(d_per_group) == 1L){
    d_per_group[[1]]
  } else {
    bdiag(d_per_group)
  }
}
