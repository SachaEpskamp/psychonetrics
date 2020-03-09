
# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_ml_lvm_group <- function(within_latent,within_residual,between_latent,between_residual,...){

  # Extract from dots things I need:
  dots <- list(...)
  design <- dots$design
  P <- dots$P
  lambda <- dots$lambda
  # lambda <- dots$lambda
  L_eta <- dots$L_eta
  # L_eta <- dots$L_eta
  L_y <- dots$L_y
  
  # Number of variables:
  nVar <- nrow(design)
  
  # Number of cases in each cluster:
  nMaxInCluster <- ncol(design)
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  
  
  # I need to construct the Jacobian for the following "observations:"
  nobs <- nVar + # Means
    (nVar * (nVar+1))/2 + # model for an individual
    nVar^2  # covariances between members of cluster
  
  # total number of elements:
  nelement <- nVar + # intercepts in nu
    nLat + # Latent intercepts
    nVar * nLat + # Factor loadings
    nLat * nLat + # Beta within
    nLat * (nLat+1) / 2 + # within factor variances
    nVar * (nVar + 1) / 2 + # Within residuals
    nLat * nLat + # Beta between
    nLat * (nLat + 1) / 2 + # Between-subject factor variances
    nVar * (nVar + 1) / 2 # between residuals
  
  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement, sparse = FALSE)
  
  # Indices:
  meanInds <- 1:nVar
  sigInds_inPerson <- nVar + seq_len(nVar*(nVar+1)/2)
  # sigInds_betPersons <-  max(sigInds_inPerson) + seq_len(nVar*(nVar+1)/2)
  sigInds_betPersons <-  max(sigInds_inPerson) + seq_len(nVar^2)
  
  # Indices model:
  nu_inds <- seq_len(nVar)
  nu_eta_inds <- max(nu_inds) + seq_len(nLat)
  lambda_inds <- max(nu_eta_inds) + seq_len(nVar * nLat)
  
  
  beta_within_inds <- max(lambda_inds) + seq_len(nLat^2)
  sigma_zeta_within_inds <- max(beta_within_inds) + seq_len(nLat * (nLat+1) / 2)
  sigma_epsilon_within_inds <- max(sigma_zeta_within_inds) + seq_len(nVar * (nVar + 1) / 2)
  
  beta_between_inds <- max(sigma_epsilon_within_inds) + seq_len(nLat^2)
  sigma_zeta_between_inds <- max(beta_between_inds) + seq_len(nLat * (nLat + 1) / 2)
  sigma_epsilon_between_inds <- max(sigma_zeta_between_inds) + seq_len(nVar * (nVar + 1) / 2)

  ### Augmentation parts (might not be needed) ###
  # Within latent
  if (within_latent == "chol"){
    aug_within_latent <-   d_sigma_cholesky(lowertri=dots$lowertri_zeta_within,L=L_eta,C=dots$C_eta_eta,In=dots$I_eta)
  } else if (within_latent == "prec"){
    aug_within_latent <- d_sigma_kappa(L = L_eta, D = dots$D_eta, sigma = dots$sigma_zeta_within)
  } else if (within_latent == "ggm"){
    aug_within_latent <- cbind(
           d_sigma_omega(L = L_eta, delta_IminOinv = dots$delta_IminOinv_zeta_within, A = dots$A_eta, delta = dots$delta_zeta_within, Dstar = dots$Dstar_eta),
           d_sigma_delta(L = L_eta,  delta_IminOinv = dots$delta_IminOinv_zeta_within,In=dots$I_eta,delta=dots$delta_zeta_within,A=dots$A_eta)
     )
  } else if (within_latent == "cov"){
    aug_within_latent <- Diagonal(nLat*(nLat+1)/2)
  }
  
  # Between latent
  if (between_latent == "chol"){
    aug_between_latent <-   d_sigma_cholesky(lowertri=dots$lowertri_zeta_between,L=L_eta,C=dots$C_eta_eta,In=dots$I_eta)
  } else if (between_latent == "prec"){
    aug_between_latent <- d_sigma_kappa(L = L_eta, D = dots$D_eta, sigma = dots$sigma_zeta_between)
  } else if (between_latent == "ggm"){
    aug_between_latent <- cbind(
      d_sigma_omega(L = L_eta, delta_IminOinv = dots$delta_IminOinv_zeta_between, A = dots$A_eta, delta = dots$delta_zeta_between, Dstar = dots$Dstar_eta),
      d_sigma_delta(L = L_eta,  delta_IminOinv = dots$delta_IminOinv_zeta_between,In=dots$I_eta,delta=dots$delta_zeta_between,A=dots$A_eta)
    )
  } else if (between_latent == "cov"){
    aug_between_latent <- Diagonal(nLat*(nLat+1)/2)
  }
  
  # Within residual
  if (within_residual == "chol"){
    aug_within_residual <-   d_sigma_cholesky(lowertri=dots$lowertri_epsilon_within,L=L_y,C=dots$C_y_y,In=dots$I_y)
  } else if (within_residual == "prec"){
    aug_within_residual <- d_sigma_kappa(L = L_y, D = dots$D_y, sigma = dots$sigma_epsilon_within)
  } else if (within_residual == "ggm"){
    aug_within_residual <- cbind(
      d_sigma_omega(L = L_y, delta_IminOinv = dots$delta_IminOinv_epsilon_within, A = dots$A_y, delta = dots$delta_epsilon_within, Dstar = dots$Dstar_y),
      d_sigma_delta(L = L_y,  delta_IminOinv = dots$delta_IminOinv_epsilon_within,In=dots$I_y,delta=dots$delta_epsilon_within,A=dots$A_y)
    )
  } else if (within_residual == "cov"){
    aug_within_residual <- Diagonal(nVar*(nVar+1)/2)
  }
  
  # Between residual
  if (between_residual == "chol"){
    aug_between_residual <-   d_sigma_cholesky(lowertri=dots$lowertri_epsilon_between,L=L_y,C=dots$C_y_y,In=dots$I_y)
  } else if (between_residual == "prec"){
    aug_between_residual <- d_sigma_kappa(L = L_y, D = dots$D_y, sigma = dots$sigma_epsilon_between)
  } else if (between_residual == "ggm"){
    aug_between_residual <- cbind(
      d_sigma_omega(L = L_y, delta_IminOinv = dots$delta_IminOinv_epsilon_between, A = dots$A_y, delta = dots$delta_epsilon_between, Dstar = dots$Dstar_y),
      d_sigma_delta(L = L_y,  delta_IminOinv = dots$delta_IminOinv_epsilon_between,In=dots$I_y,delta=dots$delta_epsilon_between,A=dots$A_y)
    )
  } else if (between_residual == "cov"){
    aug_between_residual <- Diagonal(nVar*(nVar+1)/2)
  }


  # Start filling the Jacobian:
  # fill intercept part:
  Jac[meanInds,nu_inds] <- Diagonal(nVar)
  
  # Fill latent mean part:
  Jac[meanInds,nu_eta_inds] <- d_mu_nu_eta_lvm(dots$Lambda_BetaStar_between)
  
  # Fill mean to factor loading part:
  Jac[meanInds,lambda_inds] <- d_mu_lambda_lvm(dots$nu_eta,dots$BetaStar_between,dots$I_y) 
  
  # # Fill s0 to lambda part:
  # Jac[sigInds[[1]],lambda_inds] <- d_sigmak_lambda_ml_lvm(k=1,...)
  # 
  # 
  ### Gradients for block in person ###

  # Fill s_in to lambda part:
  Jac[sigInds_inPerson,lambda_inds] <- 
    d_sigma_lambda_lvm(dots$L_y,dots$Lambda_BetaStar_within,dots$Betasta_sigmaZeta_within,dots$I_y,dots$C_y_eta,...) + 
    d_sigma_lambda_lvm(dots$L_y,dots$Lambda_BetaStar_between,dots$Betasta_sigmaZeta_between,dots$I_y,dots$C_y_eta,...)

  ### Within cluster models ###
  
  # Beta part:
  Jac[sigInds_inPerson,beta_within_inds] <- 
    d_sigma_beta_lvm(dots$L_y, dots$lambda, dots$Betasta_sigmaZeta_within, dots$C_eta_eta,dots$I_eta,dots$tBetakronBeta_within)
  
  # Fill s_in  to sigma_zeta_within part:
  Jac[sigInds_inPerson,sigma_zeta_within_inds] <- d_sigma_sigma_zeta_lvm(dots$L_y,dots$Lambda_BetaStar_within,dots$D_eta) %*% aug_within_latent
  
  # Fill s_in to sigma_epsilon_within:
  Jac[sigInds_inPerson,sigma_epsilon_within_inds] <- Diagonal(nVar * (nVar+1) / 2)  %*% aug_within_residual
  
  ### Between cluster models ###
  
  # Beta part:
  beta_between_part <- d_sigma_beta_lvm(dots$L_y, dots$lambda, dots$Betasta_sigmaZeta_between, dots$C_eta_eta,dots$I_eta,dots$tBetakronBeta_between)
  Jac[sigInds_inPerson,beta_between_inds] <- beta_between_part
    
  # Fill s_in  to sigma_zeta_between part (and store for later use):
  J_sigma_zeta_between <- d_sigma_sigma_zeta_lvm(dots$L_y,dots$Lambda_BetaStar_between,dots$D_eta) %*% aug_between_latent
  Jac[sigInds_inPerson,sigma_zeta_between_inds] <- J_sigma_zeta_between
  
  # Fill s_in to sigma_epsilon_within:
  J_sigma_epsilon_between <- Diagonal(nVar * (nVar+1) / 2)  %*% aug_between_residual
  Jac[sigInds_inPerson,sigma_epsilon_between_inds] <- J_sigma_epsilon_between
  
  ### Gradients for block between persons ###
  
  # Fill s_in to lambda part:
  Jac[sigInds_betPersons,lambda_inds] <- 
    d_sigma_lambda_lvm(Diagonal(nVar^2),dots$Lambda_BetaStar_between,dots$Betasta_sigmaZeta_between,dots$I_y,dots$C_y_eta,...)
  
    ### Between cluster models ###
  
  # Beta part:
  Jac[sigInds_betPersons,beta_between_inds] <- dots$D_y %*% beta_between_part
  
  # Fill s_in  to sigma_zeta_between part (and store for later use):
  Jac[sigInds_betPersons,sigma_zeta_between_inds] <- dots$D_y %*% J_sigma_zeta_between
  
  # Fill s_in to sigma_epsilon_within:
  Jac[sigInds_betPersons,sigma_epsilon_between_inds] <- dots$D_y %*% J_sigma_epsilon_between

  # Permute the matrix:
  Jac <- P %*% Jac
  
  # Make sparse if needed:
  # Jac <- as(Jac, "Matrix")
  Jac <- sparseordense(Jac)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_ml_lvm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_ml_lvm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




