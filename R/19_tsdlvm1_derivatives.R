d_mu_mu_eta_tsdlvm1 <- function(lambda,...){
  lambda
}

d_mu_lambda_tsdlvm1 <- function(mu_eta,I_y,...){
  t(mu_eta) %x% I_y
}

d_sigmak_lambda_tsdlvm1 <- function(lambda, k = 0, Sigma_eta_0, Sigma_eta_1, C_y_eta, I_y, L_y, ...){
  
  if (k == 0){
    sigEta <- Sigma_eta_0
  } else {
    sigEta <- Sigma_eta_1
  }
  
  within <- (
    (I_y %x% (lambda %*% sigEta)) %*% C_y_eta + 
      ( (lambda %*% t(sigEta)) %x%  I_y)  
  )
  
  res <- within 
  
  if (k == 0){
    return(L_y %*% res)
  } else {
    return(res)
  }
}

d_sigma0_sigma_zeta_tsdlvm1 <- function(lamWkronlamW, BetaStar, D_eta, ...){
  # lamWkronlamW %*% BetaStar %*% D_eta
  BetaStar %*% D_eta
}

d_sigma0_beta_tsdlvm1 <- function(BetaStar,I_eta,Sigma_eta_1,C_eta_eta,...){
  
  res <- ((I_eta%x%I_eta) + C_eta_eta) %*% BetaStar %*% (Sigma_eta_1 %x% I_eta) 
 
  return(res)

}

d_sigma1_beta_tsdlvm1 <- function(J_sigma_beta,IkronBeta,k,  Sigma_eta_0,I_eta,...){
  
  IkronBeta %*% J_sigma_beta + (t(Sigma_eta_0) %x% I_eta)
  
}


# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_tsdlvm1_group <- function(zeta,epsilon,...){
  # Extract from dots things I need:
  dots <- list(...)
  P <- dots$P
  lambda <- dots$lambda
  
  # lambda <- dots$lambda
  L_eta <- dots$L_eta
  # L_eta <- dots$L_eta
  L_y <- dots$L_y
  
  
  # Number of variables:
  nvar <- nrow(lambda) * 2
  
  # Number of nodes:
  nNode <- nrow(lambda)
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nNode + # Exogenous means
    nNode + # intercepts
    nLat + # Latent means
    nNode*(nNode+1)/2 + # Exogenous variances
    nNode * nLat + # Factor loadings
    nLat * (nLat + 1) / 2 + # Contemporaneous
    nLat^2 + # Beta
    nNode * (nNode+1) / 2 # Residual
  
  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement, sparse = FALSE)
  
  # Indices:
  meanInds_exo <- 1:nNode
  meanInds_endo <- nNode + 1:nNode
  sigmaStarInds <- nvar + seq_len(nNode*(nNode+1)/2)
  sigma0Inds <- max(sigmaStarInds) + seq_len(nNode*(nNode+1)/2)
  sigma1Inds <- max(sigma0Inds) + seq_len(nNode^2)
  
  # Indices model:
  exomeanInds <- 1:nNode
  interceptInds <- max(exomeanInds) + 1:nNode
  latmeanInds <- max(interceptInds) + 1:nLat
  exovarInds <-  max(latmeanInds) + seq_len(nNode*(nNode+1)/2)
  lambdaInds <- max(exovarInds) + seq_len(nNode * nLat)
  contInds <- max(lambdaInds) + seq_len(nLat * (nLat + 1) / 2)
  betaInds <- max(contInds) + seq_len(nLat^2)
  residInds <- max(betaInds) +  seq_len(nNode*(nNode+1)/2)
  
  ### Augmentation parts ###
  
  # Contemporaneous
  if (zeta == "chol"){
    aug_zeta <-   d_sigma_cholesky(lowertri=dots$lowertri_zeta,L=L_eta,C=dots$C_eta_eta,In=dots$I_eta)
  } else if (zeta == "prec"){
    aug_zeta <- d_sigma_kappa(L = L_eta, D = dots$D_eta, sigma = dots$sigma_zeta)
  } else if (zeta == "ggm"){
    aug_zeta <- cbind(
      d_sigma_omega(L = L_eta, delta_IminOinv = dots$delta_IminOinv_zeta, A = dots$A_eta, delta = dots$delta_zeta, Dstar = dots$Dstar_eta),
      d_sigma_delta(L = L_eta,  delta_IminOinv = dots$delta_IminOinv_zeta,In=dots$I_eta,delta=dots$delta_zeta,A=dots$A_eta)
    )
  } else if (zeta == "cov"){
    aug_zeta <- Diagonal(nLat*(nLat+1)/2)
  }
  
  # Residual
  if (epsilon == "chol"){
    aug_epsilon <-   d_sigma_cholesky(lowertri=dots$lowertri_epsilon,L=L_y,C=dots$C_y_y,In=dots$I_y)
  } else if (epsilon == "prec"){
    aug_epsilon <- d_sigma_kappa(L = L_y, D = dots$D_y, sigma = dots$sigma_epsilon)
  } else if (epsilon == "ggm"){
    aug_epsilon <- cbind(
      d_sigma_omega(L = L_y, delta_IminOinv = dots$delta_IminOinv_epsilon, A = dots$A_y, delta = dots$delta_epsilon, Dstar = dots$Dstar_y),
      d_sigma_delta(L = L_y,  delta_IminOinv = dots$delta_IminOinv_epsilon,In=dots$I_y,delta=dots$delta_epsilon,A=dots$A_y)
    )
  } else if (epsilon == "cov"){
    aug_epsilon <- Diagonal(nNode*(nNode+1)/2)
  }
  
  # Exogenous mean part:
  Jac[meanInds_exo, exomeanInds] <- Diagonal(nNode)
  
  
  # Intercept part:
  Jac[meanInds_endo, interceptInds] <- Diagonal(nNode)
  
  # Fill latent mean part:
  Jac[meanInds_endo,latmeanInds] <- d_mu_mu_eta_tsdlvm1( ...)
  
  # Fill mean to factor loading part:
  Jac[meanInds_endo,lambdaInds] <- d_mu_lambda_tsdlvm1(...)
  
  
  
  # Exogenous block variances:
  Jac[sigmaStarInds,exovarInds] <- d_sigma_cholesky(lowertri=dots$exo_cholesky,L=L_y,C=dots$C_y_y,In=dots$I_y)
  
  ### Sigma 0 part:
  
  # Sigma 0 to lambda:
  Jac[sigma0Inds,lambdaInds] <- d_sigmak_lambda_tsdlvm1(k=0,...)
  
  # Sigma 0 to contemporaneous
  J_sigma_zeta <- d_sigma0_sigma_zeta_tsdlvm1(...) %*% aug_zeta
  Jac[sigma0Inds,contInds] <- L_y %*% dots$lamWkronlamW %*% J_sigma_zeta 
  
  # Fill s0 to beta part (and store for later use):
  J_sigma_beta <- d_sigma0_beta_tsdlvm1(...)
  Jac[sigma0Inds, betaInds] <- L_y  %*% dots$lamWkronlamW %*% J_sigma_beta
  
  # Fill s0 to sigma_epsilon_within:
  Jac[sigma0Inds, residInds] <- Diagonal(nNode * (nNode+1) / 2)  %*% aug_epsilon
  
  
  ## Sigma 1 parts
  
  # Fill s1 to lambda part:
  Jac[sigma1Inds, lambdaInds] <- d_sigmak_lambda_tsdlvm1(k=1,...)
  
  
  
  # Fill s1 to sigma_zeta part:
  J_sigma_zeta <- dots$IkronBeta %*% J_sigma_zeta
  Jac[sigma1Inds, contInds] <- dots$lamWkronlamW %*% J_sigma_zeta
  
  # Fill s1 to beta part (and store for later use):
  J_sigma_beta <- d_sigma1_beta_tsdlvm1(J_sigma_beta=J_sigma_beta,...)
  Jac[sigma1Inds,betaInds] <- dots$lamWkronlamW %*% J_sigma_beta
  
  
  
  # Permute the matrix:
  Jac <- P %*% Jac

  # Make sparse if needed:
  Jac <- as(Jac, "Matrix")
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_tsdlvm1 <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_tsdlvm1_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




