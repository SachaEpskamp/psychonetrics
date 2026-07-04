# Derivatives for the panelvar model. These are the dlvm1 derivatives with
# lambda = I: no factor loadings, no residual (co)variances, and observed
# means in place of latent means. See 18_dlvm1_derivatives.R for the general
# versions.

# Without elimination matrix (maps vech(sigma_zeta_within) to vec(sigma_0)):
d_sigma0_sigma_zeta_within_panelvar <- function(BetaStar, D2, ...){
  BetaStar %*% D2
}

# Without elimination matrix (maps vec(beta) to vec(sigma_0)):
d_sigma0_beta_panelvar <- function(BetaStar,In,allSigmas_within,C,...){
  ((In %x% In) + C) %*% BetaStar %*% (allSigmas_within[[2]] %x% In)
}

# Lag-k (k > 1) part (maps vec(beta) to vec(sigma_k)):
d_sigmak_beta_panelvar <- function(J_sigma_beta,IkronBeta,k,allSigmas_within,In,...){
  IkronBeta %*% J_sigma_beta + (t(allSigmas_within[[k-1]]) %x% In)
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_panelvar_group <- function(within_latent,between_latent,...){

  # Extract from dots things I need:
  dots <- list(...)
  design <- dots$design
  P <- dots$P
  L <- dots$L
  In <- dots$In

  # Number of variables:
  nVar <- nrow(design)

  # Number of time points:
  nTime <- ncol(design)

  # I need to construct the Jacobian for the following "observations:"
  nobs <- nVar + # Means
    (nVar * (nVar+1))/2 + # Variances
    nVar^2 * (nTime - 1) # lagged variances

  # total number of elements:
  nelement <- nVar + # means
    nVar * (nVar+1) / 2 + # within-person contemporaneous block
    nVar^2 + # temporal effects
    nVar * (nVar + 1) / 2 # between-person block

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement, sparse = FALSE)

  # Indices (observations):
  meanInds <- 1:nVar
  sigInds <- list(
    nVar + seq_len(nVar*(nVar+1)/2)
  )

  # For each lag:
  for (t in 2:nTime){
    sigInds[[t]] <- max(sigInds[[t-1]]) + seq_len(nVar^2)
  }

  # Indices model:
  mu_inds <- seq_len(nVar)
  sigma_zeta_within_inds <- max(mu_inds) + seq_len(nVar * (nVar+1) / 2)
  beta_inds <- max(sigma_zeta_within_inds) + seq_len(nVar^2)
  sigma_zeta_between_inds <- max(beta_inds) + seq_len(nVar * (nVar + 1) / 2)

  ### Augmentation parts ###

  # Within-person block:
  if (within_latent == "chol"){
    aug_within_latent <- d_sigma_cholesky(lowertri=dots$lowertri_zeta_within,L=L,C=dots$C,In=In)
  } else if (within_latent == "prec"){
    aug_within_latent <- d_sigma_kappa(L = L, D = dots$D2, sigma = dots$sigma_zeta_within)
  } else if (within_latent == "ggm"){
    aug_within_latent <- cbind(
      d_sigma_omega(L = L, delta_IminOinv = dots$delta_IminOinv_zeta_within, A = dots$A, delta = dots$delta_zeta_within, Dstar = dots$Dstar),
      d_sigma_delta(L = L,  delta_IminOinv = dots$delta_IminOinv_zeta_within,In=In,delta=dots$delta_zeta_within,A=dots$A)
    )
  } else if (within_latent == "cov"){
    aug_within_latent <- Diagonal(nVar*(nVar+1)/2)
  } else if (within_latent == "cor"){
    aug_within_latent <- cbind(
      d_sigma_rho(L = L, SD = dots$SD_zeta_within, A = dots$A, Dstar = dots$Dstar),
      d_sigma_SD(L = L, SD_IplusRho = dots$SD_IplusRho_zeta_within, In = In, A = dots$A)
    )
  }

  # Between-person block:
  if (between_latent == "chol"){
    aug_between_latent <- d_sigma_cholesky(lowertri=dots$lowertri_zeta_between,L=L,C=dots$C,In=In)
  } else if (between_latent == "prec"){
    aug_between_latent <- d_sigma_kappa(L = L, D = dots$D2, sigma = dots$sigma_zeta_between)
  } else if (between_latent == "ggm"){
    aug_between_latent <- cbind(
      d_sigma_omega(L = L, delta_IminOinv = dots$delta_IminOinv_zeta_between, A = dots$A, delta = dots$delta_zeta_between, Dstar = dots$Dstar),
      d_sigma_delta(L = L,  delta_IminOinv = dots$delta_IminOinv_zeta_between,In=In,delta=dots$delta_zeta_between,A=dots$A)
    )
  } else if (between_latent == "cov"){
    aug_between_latent <- Diagonal(nVar*(nVar+1)/2)
  } else if (between_latent == "cor"){
    aug_between_latent <- cbind(
      d_sigma_rho(L = L, SD = dots$SD_zeta_between, A = dots$A, Dstar = dots$Dstar),
      d_sigma_SD(L = L, SD_IplusRho = dots$SD_IplusRho_zeta_between, In = In, A = dots$A)
    )
  }

  # Fill mean part:
  Jac[meanInds,mu_inds] <- Diagonal(nVar)

  # Fill s0 to sigma_zeta_within part (and store for later use):
  J_sigma_zeta_within <- d_sigma0_sigma_zeta_within_panelvar(...) %*% aug_within_latent
  Jac[sigInds[[1]],sigma_zeta_within_inds] <- L %*% J_sigma_zeta_within

  # Fill s0 to beta part (and store for later use):
  J_sigma_beta <- d_sigma0_beta_panelvar(...)
  Jac[sigInds[[1]],beta_inds] <- L %*% J_sigma_beta

  # Fill s0 to sigma_zeta_between part (and store for later use). Note that
  # L %*% D2 = I, so the s0 rows are the augmentation matrix itself:
  J_sigma_zeta_between <- dots$D2 %*% aug_between_latent
  Jac[sigInds[[1]],sigma_zeta_between_inds] <- L %*% J_sigma_zeta_between

  # For every further lag:
  for (t in 2:nTime){
    # Fill sk to sigma_zeta_within part (and store for later use):
    J_sigma_zeta_within <- dots$IkronBeta %*% J_sigma_zeta_within
    Jac[sigInds[[t]],sigma_zeta_within_inds] <- J_sigma_zeta_within

    # Fill sk to beta part (and store for later use):
    J_sigma_beta <- d_sigmak_beta_panelvar(J_sigma_beta=J_sigma_beta,k=t,...)
    Jac[sigInds[[t]],beta_inds] <- J_sigma_beta

    # Fill sk to sigma_zeta_between part:
    Jac[sigInds[[t]],sigma_zeta_between_inds] <- J_sigma_zeta_between
  }

  # Permute the matrix:
  Jac <- P %*% Jac

  # Make sparse if needed:
  Jac <- as(Jac, "Matrix")

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_panelvar <- function(prep){
  # model is already prepared!

  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_panelvar_group)

  # Bind by column and return:
  Reduce("bdiag",d_per_group)
}
