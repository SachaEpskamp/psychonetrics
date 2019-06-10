d_mu_mu_eta_dlvm1 <- function(lambda,...){
  lambda
}

# d_mu_lambda_dlvm1 <- function(mu_eta,I_y,...){
#   t(mu_eta) %x% I_y
# }

d_mu_lambda_dlvm1 <- function(mu_eta,I_y,...){
  t(mu_eta) %x% I_y
}

# d_sigmak_lambda_dlvm1 <- function(lambda, k = 0,  allSigmas_within, C_y_eta, I_y, L_y, ...){
#   res <- (
#     (I_y %x% (lambda %*% allSigmas_within[[k]])) %*% C_y_eta + 
#     ( (lambda %*% t(allSigmas_within[[k]])) %x%  I_y)  
#   )
#   # browser()
#   
#   
#   # res2 <- ((I_y %x% I_y) + C_y_eta) %*% ((lambda %*% allSigmas_within[[k]]) %x%  I_y) 
#   
#   if (k == 1){
#     return(L_y %*% res)
#   } else {
#     return(res)
#   }
# }

d_sigmak_lambda_dlvm1 <- function(lambda, k = 0,  allSigmas_within, C_y_eta, I_y, L_y,sigma_zeta_between, ...){
  within <- (
    (I_y %x% (lambda %*% allSigmas_within[[k]])) %*% C_y_eta + 
      ( (lambda %*% t(allSigmas_within[[k]])) %x%  I_y)  
  )
  
  between <- (I_y %x% (lambda %*% sigma_zeta_between)) %*% C_y_eta + 
    ( (lambda %*% sigma_zeta_between) %x%  I_y)  
  
  res <- within + between
  
  if (k == 1){
    return(L_y %*% res)
  } else {
    return(res)
  }
}
# 
# d_sigmak_lambda_dlvm1 <- function(lambda, sigma_zeta_between, C_y_between, I_y, ...){
#   (I_y %x% (lambda %*% sigma_zeta_between)) %*% C_y_between + 
#     ( (lambda %*% sigma_zeta_between) %x%  I_y)  
# }

# Without elimination matrix
d_sigma0_sigma_zeta_within_dlvm1 <- function(lamWkronlamW, BetaStar, D_eta, ...){
  # lamWkronlamW %*% BetaStar %*% D_eta
  BetaStar %*% D_eta
}

# Without elimination matrix:
d_sigma0_beta_dlvm1 <- function(BetaStar,I_eta,allSigmas_within,C_eta_eta,...){
  # lamWkronlamW %*% (t(Vec(sigma_zeta_within)) %x% (I_eta %x% I_eta)) %*% (t(BetaStar) %x% BetaStar) %*% (
  #   E %x% I_eta + I_eta %x% E
  # )
  # Number of columns in betastar:
  # 
  # cp <- tcrossprod(t(Vec(sigma_zeta_within)), BetaStar)
  # res <- (cp %x% BetaStar)  %*% (
  #   E %x% I_eta + I_eta %x% E
  # )
  # 
  # res2 <- BetaStar %*% (allSigmas_within[[2]] %x% I_eta) +
  #   BetaStar %*% (I_eta %x% allSigmas_within[[2]]) %*% C_eta_eta
  # 

  res <- ((I_eta%x%I_eta) + C_eta_eta) %*% BetaStar %*% (allSigmas_within[[2]] %x% I_eta) 

    
  
  return(res)
  # # Old experimental stuff:
  # nc <- ncol(BetaStar)
  # # If the number is not *too* big, do normal:
  # if (nc < 12){
  #   res <- (t(Vec(sigma_zeta_within)) %x% (I_eta %x% I_eta)) %*% (t(BetaStar) %x% BetaStar) %*% (
  #     E %x% I_eta + I_eta %x% E
  #   )
  #   
  #   return(res)
  # } else {
  #   # Fancy method:
  #   smat <- (
  #     E %x% I_eta + I_eta %x% E
  #   )
  #   postmat <- do.call(cbind,lapply(1:nc,function(i)Vec(BetaStar%*%Matrix(smat[,i],nc,nc)%*%BetaStar)))
  #   res <- (t(Vec(sigma_zeta_within)) %x% (I_eta %x% I_eta)) %*% postmat
  #   return(res)
  # }
}

d_sigmak_beta_dlvm1 <- function(J_sigma_beta,IkronBeta,k,  allSigmas_within,I_eta,...){
  
  IkronBeta %*% J_sigma_beta + (t(allSigmas_within[[k-1]]) %x% I_eta)
  
}


# d_sigmak_lambda_dlvm1 <- function(lambda, sigma_zeta_between, C_y_eta, I_y, ...){
#     (I_y %x% (lambda %*% sigma_zeta_between)) %*% C_y_eta + 
#       ( (lambda %*% sigma_zeta_between) %x%  I_y)  
# }

d_sigmak_sigma_zeta_between_dlvm1 <- function(lambda,D_eta,...){
  (lambda %x% lambda) %*% D_eta
}

# 
#  # derivative of sigma0 (full vec) with respect to beta:
# d_sigma0_beta_dlvm1 <- function(BetaStar,E,In,w,B_Sigma_mu,C,...){
#   (t(w) %x% (In %x% In)) %*% 
#     (t(BetaStar) %x% BetaStar) %*% 
#     (E %x% In + In %x% E) - 
#     BetaStar %*% (
#       (In %x% B_Sigma_mu) %*% C + (B_Sigma_mu %x% In)
#     )
# }
# # 
# # # derivative of sigma0 with respect to sigma_zeta:
# d_sigma0_sigma_zeta_dlvm1 <- function(BetaStar,D2,...){
#   BetaStar %*% D2
# }
# 
# ## Lag k derivatives:
# d_sigmak_beta_dlvm1 <- function(k,IkronBeta,Jb,allSigmas,Sigma_mu_kron_I,In,...){
#   IkronBeta %*% Jb + (t(allSigmas[[k-1]]) %x% In) - Sigma_mu_kron_I
# }
# 
# d_sigmak_sigma_zeta_dlvm1 <- function(IkronBeta,Jsigz,...){
#   IkronBeta %*% Jsigz
# }
# 
# d_sigmak_sigma_mu_dlvm1 <- function(D2,...){
#   D2
# }
# 
### Augmentation parts:

# Derivative of sigma_zeta to cholesky:
# d_sigma_zeta_within_cholesky_dlvm1 <- function(lowertri_zeta_within,L,C_eta_eta,I_eta,...){
#   d_sigma_cholesky(lowertri=lowertri_zeta_within,L=L,C=C_eta_eta,In=I_eta)
# }
# 
# # Derivative of sigma_zeta to precision:
# d_sigma_zeta_kappa_dlvm1 <- function(L,D2, sigma_zeta,...){
#   d_sigma_kappa(L = L, D = D2, sigma = sigma_zeta)
# }
# 
# # Derivative of sigma_zeta to ggm:
# d_sigma_zeta_ggm_dlvm1 <- function(L,delta_IminOinv_zeta,A,delta_zeta,Dstar,In,...){
#   cbind(
#     d_sigma_omega(L = L, delta_IminOinv = delta_IminOinv_zeta, A = A, delta = delta_zeta, Dstar = Dstar),
#     d_sigma_delta(L = L,  delta_IminOinv = delta_IminOinv_zeta,In=In,delta=delta_zeta,A=A)
#   )
# }
# 
# 
# d_sigma_mu_cholesky_dlvm1 <- function(lowertri_mu,L,C,In,...){
#   d_sigma_cholesky(lowertri=lowertri_mu,L=L,C=C,In=In)
# }
# 
# # Derivative of sigma_mu to precision:
# d_sigma_mu_kappa_dlvm1 <- function(L,D2, sigma_mu,...){
#   d_sigma_kappa(L = L, D = D2, sigma = sigma_mu)
# }
# 
# # Derivative of sigma_mu to ggm:
# d_sigma_mu_ggm_dlvm1 <- function(L,delta_IminOinv_mu,A,delta_mu,Dstar,In,...){
#   cbind(
#     d_sigma_omega(L = L, delta_IminOinv = delta_IminOinv_mu, A = A, delta = delta_mu, Dstar = Dstar),
#     d_sigma_delta(L = L,  delta_IminOinv = delta_IminOinv_mu,In=In,delta=delta_mu,A=A)
#   )
# }

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_dlvm1_group <- function(within_latent,within_residual,between_latent,between_residual,...){
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
  
  # Number of times points:
  nTime <- ncol(design)
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  # Number of within latents:
  # nLat <- ncol(lambda)
  
  # Number of between latents:
  # nLat <- ncol(lambda)
  
  # I need to construct the Jacobian for the following "observations:"
  nobs <- nVar + # Means
    (nVar * (nVar+1))/2 + # Variances
    nVar^2 * (nTime - 1) # lagged variances
  
  # total number of elements:
  nelement <- nVar + # intercepts in tau
    nLat + # Latent means
    nVar * nLat + # Factor loadings
    nLat * (nLat+1) / 2 + # within factor variances
    nLat^2 + # temporal effects
    nVar * (nVar + 1) / 2 + # Within residuals
    # nVar * nLat + # Between-subject factor loadings
    nLat * (nLat + 1) / 2 + # Between-subject factor variances
    nVar * (nVar + 1) / 2 # between residuals
  
  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement, sparse = FALSE)
  
  # Indices:
  meanInds <- 1:nVar
  sigInds <- list(
    nVar + seq_len(nVar*(nVar+1)/2)
  )
  
  # For each lag:
  for (t in 2:nTime){
    sigInds[[t]] <- max(sigInds[[t-1]]) + seq_len(nVar^2)
  }
  
  # Indices model:
  tau_inds <- seq_len(nVar)
  mu_eta_inds <- max(tau_inds) + seq_len(nLat)
  lambda_inds <- max(mu_eta_inds) + seq_len(nVar * nLat)
  
  # lambda_inds <- max(mu_eta_inds) + seq_len(nVar * nLat)
  sigma_zeta_within_inds <- max(lambda_inds) + seq_len(nLat * (nLat+1) / 2)
  beta_inds <- max(sigma_zeta_within_inds) + seq_len(nLat^2)
  sigma_epsilon_within_inds <- max(beta_inds) + seq_len(nVar * (nVar + 1) / 2)
  
  sigma_zeta_between_inds <- max(sigma_epsilon_within_inds) + seq_len(nLat * (nLat + 1) / 2)
  sigma_epsilon_between_inds <- max(sigma_zeta_between_inds) + seq_len(nVar * (nVar + 1) / 2)
  
  ### Augmentation parts ###

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

 
  # message("Starting...")
  # fill intercept part:
  Jac[meanInds,tau_inds] <- Diagonal(nVar)
  
  # Fill latent mean part:
  Jac[meanInds,mu_eta_inds] <- d_mu_mu_eta_dlvm1( ...)
  
  # Fill mean to factor loading part:
  Jac[meanInds,lambda_inds] <- d_mu_lambda_dlvm1(...)
  
  # # Fill s0 to lambda part:
  # Jac[sigInds[[1]],lambda_inds] <- d_sigmak_lambda_dlvm1(k=1,...)
  # 
  # 
  # Fill s0 to lambda part:
  Jac[sigInds[[1]],lambda_inds] <- d_sigmak_lambda_dlvm1(k=1,...)
  
  
  # Fill s0  to sigma_zeta_within part (and store for later use):
  J_sigma_zeta_within <- d_sigma0_sigma_zeta_within_dlvm1(...) %*% aug_within_latent
  Jac[sigInds[[1]],sigma_zeta_within_inds] <- L_y %*% dots$lamWkronlamW %*% J_sigma_zeta_within 
  
  
  # Fill s0 to beta part (and store for later use):
  J_sigma_beta <- d_sigma0_beta_dlvm1(...)
  Jac[sigInds[[1]],beta_inds] <- L_y  %*% dots$lamWkronlamW %*% J_sigma_beta
  
  # Fill s0 to sigma_epsilon_within:
  Jac[sigInds[[1]],sigma_epsilon_within_inds] <- Diagonal(nVar * (nVar+1) / 2)  %*% aug_within_residual
  
  
  # Fill s0 to lambda part and store for later use
  # J_sigmak_lambda <- d_sigmak_lambda_dlvm1(...)
  # Jac[sigInds[[1]],lambda_inds] <- L_y %*% J_sigmak_lambda
  
  # Fill s0 to sigma_zeta_between, and store for later use:
  J_sigmak_sigma_zeta_between <- d_sigmak_sigma_zeta_between_dlvm1(...)  %*% aug_between_latent
  Jac[sigInds[[1]],sigma_zeta_between_inds] <- L_y %*% J_sigmak_sigma_zeta_between
   
  # Fill sigma_epsilon_between inds:
  Jac[sigInds[[1]],sigma_epsilon_between_inds] <-  aug_between_residual

  # For every further lag:
  for (t in 2:nTime){
    # message(paste0("t = ",t))
    # Fill sk to lambda part:
    Jac[sigInds[[t]],lambda_inds] <- d_sigmak_lambda_dlvm1(k=t,...)

    # Fill sk  to sigma_zeta_within part (and store for later use):
    J_sigma_zeta_within <- dots$IkronBeta %*% J_sigma_zeta_within
    Jac[sigInds[[t]],sigma_zeta_within_inds] <- dots$lamWkronlamW %*% J_sigma_zeta_within
    
    # Fill sk to beta part (and store for later use):
    J_sigma_beta <- d_sigmak_beta_dlvm1(J_sigma_beta=J_sigma_beta,k=t,...)
    Jac[sigInds[[t]],beta_inds] <- dots$lamWkronlamW %*% J_sigma_beta
    # 
    # # Fill sk to sigma_epsilon_within:
    # Jac[sigInds[[1]],sigma_epsilon_within_inds] <- Diagonal(nVar * (nVar+1) / 2)
    # 
    
    # Fill sk to lambda part and store for later use
    # Jac[sigInds[[t]],lambda_inds] <- J_sigmak_lambda
    
    # Fill s0 to sigma_zeta_between, and store for later use:
    Jac[sigInds[[t]],sigma_zeta_between_inds] <- J_sigmak_sigma_zeta_between
    
    # Fill sigma_epsilon_between inds:
    Jac[sigInds[[t]],sigma_epsilon_between_inds] <- dots$D_y %*% aug_between_residual
  }
  
  # Permute the matrix:
  Jac <- P %*% Jac
  
  # Make sparse if needed:
  Jac <- as(Jac, "Matrix")

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_dlvm1 <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_dlvm1_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




