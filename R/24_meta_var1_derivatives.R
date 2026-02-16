# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_meta_var1_group <- function(beta, zeta, randomEffects, cpp, ...){

  dots <- list(...)

  # Number of nodes:
  nNode <- nrow(beta)

  # Number of modeled elements in the mean part:
  # vech(Sigma0) = nNode*(nNode+1)/2 + vec(Sigma1) = nNode^2
  nSigma0 <- nNode * (nNode + 1) / 2
  nSigma1 <- nNode^2
  nmod <- nSigma0 + nSigma1

  # Number of observations: mean part + variance part
  nobs <- nmod +                  # mean part (vech(Sigma0), vec(Sigma1))
    nmod * (nmod + 1) / 2         # variance part (vech of random effects covariance)

  # Mean part and variance part indices:
  meanPart <- seq_len(nmod)
  varPart <- max(meanPart) + seq_len(nmod * (nmod + 1) / 2)

  # Sigma0 and Sigma1 row indices within mean part:
  sigma0Inds <- seq_len(nSigma0)
  sigma1Inds <- nSigma0 + seq_len(nSigma1)

  #### VAR(1) parameter count ####
  # Parameters: beta (nNode^2) + sigma_zeta (nNode*(nNode+1)/2)
  nBeta <- nNode^2
  nSigmaZeta <- nNode * (nNode + 1) / 2

  # Adjust sigma_zeta parameter count for ggm/cor parameterizations:
  if (zeta == "ggm"){
    nSigmaZeta <- nNode * (nNode - 1) / 2 + nNode  # omega (off-diag) + delta (diag)
  } else if (zeta == "prec"){
    nSigmaZeta <- nNode * (nNode + 1) / 2
  }
  nVAR1 <- nBeta + nSigmaZeta

  #### Random effects parameter count ####
  if (randomEffects == "cov"){
    nRan <- nmod * (nmod + 1) / 2
  } else if (randomEffects == "chol"){
    nRan <- nmod * (nmod + 1) / 2
  } else if (randomEffects == "prec"){
    nRan <- nmod * (nmod + 1) / 2
  } else if (randomEffects == "ggm"){
    nRan <- nmod * (nmod - 1) / 2 + nmod
  } else if (randomEffects == "cor"){
    nRan <- nmod * (nmod - 1) / 2 + nmod
  }

  # Total parameters:
  nTotal <- nVAR1 + nRan

  # Empty Jacobian:
  Jac <- matrix(0, nobs, nTotal)

  # Parameter column indices:
  betaInds_col <- seq_len(nBeta)
  sigmazetaInds_col <- nBeta + seq_len(nSigmaZeta)
  ranInds_col <- nVAR1 + seq_len(nRan)

  #### Fill mean part: VAR(1) Jacobian ####
  # We need a dummy full sigma (2n x 2n) for the _cpp helper functions
  # that extract sigma1 = sigma[n+(1:n), 1:n] and sigma0 = sigma[n+(1:n), n+(1:n)]
  dummySigma <- matrix(0, 2*nNode, 2*nNode)
  dummySigma[nNode + (1:nNode), nNode + (1:nNode)] <- as.matrix(dots$Sigma0)
  dummySigma[nNode + (1:nNode), 1:nNode] <- as.matrix(dots$Sigma1)
  dummySigma[1:nNode, nNode + (1:nNode)] <- t(as.matrix(dots$Sigma1))

  # d(vech(Sigma0))/d(beta):
  Jac[sigma0Inds, betaInds_col] <- Jb <- d_sigma0_beta_var1_cpp(
    BetaStar = dots$BetaStar, In = dots$In, sigma = dummySigma, C = dots$C, L = dots$L
  )

  # d(vech(Sigma0))/d(sigma_zeta):
  Jac[sigma0Inds, sigmazetaInds_col] <- d_sigma0_sigma_zeta_var1_cpp(
    L = dots$L, BetaStar = dots$BetaStar, D2 = dots$D2
  )

  # Chain rule for sigma_zeta parameterization:
  if (zeta == "chol"){
    Jac[sigma0Inds, sigmazetaInds_col] <- Jac[sigma0Inds, sigmazetaInds_col] %*%
      d_sigma_zeta_cholesky_var1_cpp(
        lowertri_zeta = dots$lowertri_zeta, L = dots$L, C = dots$C, In = dots$In
      )
  } else if (zeta == "prec"){
    Jac[sigma0Inds, sigmazetaInds_col] <- Jac[sigma0Inds, sigmazetaInds_col] %*%
      d_sigma_zeta_kappa_var1_cpp(
        L = dots$L, D2 = dots$D2, sigma_zeta = dots$sigma_zeta
      )
  } else if (zeta == "ggm"){
    Jac[sigma0Inds, sigmazetaInds_col] <- Jac[sigma0Inds, sigmazetaInds_col] %*%
      d_sigma_zeta_ggm_var1_cpp(
        L = dots$L, delta_IminOinv_zeta = dots$delta_IminOinv_zeta, A = dots$A,
        delta_zeta = dots$delta_zeta, Dstar = dots$Dstar, In = dots$In
      )
  }

  # Store for sigma1 derivatives:
  Js <- Jac[sigma0Inds, sigmazetaInds_col]

  # d(vec(Sigma1))/d(beta):
  Jac[sigma1Inds, betaInds_col] <- d_sigma1_beta_var1_cpp(
    IkronBeta = dots$IkronBeta, D2 = dots$D2, Jb = Jb,
    sigma = dummySigma, beta = beta, In = dots$In
  )

  # d(vec(Sigma1))/d(sigma_zeta):
  Jac[sigma1Inds, sigmazetaInds_col] <- d_sigma1_sigma_zeta_var1_cpp(
    IkronBeta = dots$IkronBeta, D2 = dots$D2, Js = Js
  )

  #### Fill variance part: random effects Jacobian ####
  # Same structure as meta_lvm
  if (randomEffects == "cov"){
    Jac[varPart, ranInds_col] <- as.matrix(Diagonal(nRan))

  } else if (randomEffects == "chol"){
    if (cpp){
      Jac[varPart, ranInds_col] <- d_sigma_cholesky_cpp(
        lowertri = dots$lowertri_randomEffects, L = dots$L_c, C = dots$C_c, In = dots$In_c
      )
    } else {
      Jac[varPart, ranInds_col] <- d_sigma_cholesky(
        lowertri = dots$lowertri_randomEffects, L = dots$L_c, C = dots$C_c, In = dots$In_c
      )
    }

  } else if (randomEffects == "ggm"){
    netPart <- nVAR1 + seq_len(nmod * (nmod - 1) / 2)
    scalingPart <- max(netPart) + seq_len(nmod)

    if (cpp){
      Jac[varPart, netPart] <- d_sigma_omega_cpp(
        L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c
      )
      Jac[varPart, scalingPart] <- d_sigma_delta_cpp(
        L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        In = dots$In_c, A = dots$A_c
      )
    } else {
      Jac[varPart, netPart] <- d_sigma_omega(
        L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c
      )
      Jac[varPart, scalingPart] <- d_sigma_delta(
        L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        In = dots$In_c, A = dots$A_c
      )
    }

  } else if (randomEffects == "prec"){
    if (cpp){
      Jac[varPart, ranInds_col] <- d_sigma_kappa_cpp(
        L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects
      )
    } else {
      Jac[varPart, ranInds_col] <- d_sigma_kappa(
        L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects
      )
    }

  } else if (randomEffects == "cor"){
    corPart <- nVAR1 + seq_len(nmod * (nmod - 1) / 2)
    sdPart <- max(corPart) + seq_len(nmod)

    if (cpp){
      Jac[varPart, corPart] <- d_sigma_rho_cpp(
        L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c
      )
      Jac[varPart, sdPart] <- d_sigma_SD_cpp(
        L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects,
        In = dots$In_c, A = dots$A_c
      )
    } else {
      Jac[varPart, corPart] <- d_sigma_rho(
        L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c
      )
      Jac[varPart, sdPart] <- d_sigma_SD(
        L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects,
        In = dots$In_c, A = dots$A_c
      )
    }
  }

  # Make sparse if needed:
  Jac <- sparseordense(Jac)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_meta_var1 <- function(prep){
  # model is already prepared!

  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels, do.call, what = d_phi_theta_meta_var1_group)

  # Bind by column and return:
  Reduce("bdiag", d_per_group)
}
