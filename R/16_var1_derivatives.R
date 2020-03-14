# Derivative of nu with respect to mu:
d_mu_mu_var1 <- function(beta,...){
  Diagonal(nrow(beta))
}

# Derivative of exogenous variances part
d_sigmastar_exo_cholesky_var1 <- function(In, L, C, exo_cholesky, ...){
  L %*% (
    ( In %x% exo_cholesky) %*% C %*% t(L) + 
      (exo_cholesky %x% In) %*% t(L)
  )
}

# derivative of sigma0 with respect to beta:
d_sigma0_beta_var1 <- function(BetaStar,In, sigma,C, L,...){

  n <- nrow(In)
  sigma1 <- sigma[n + (1:n), (1:n)]
  res <- L %*% ((In%x%In) + C) %*% BetaStar %*% (sigma1 %x% In) 
  return(res)
  
  # ((t(sigmaZetaVec) %*% t(BetaStar)) %x% L_betaStar) %*% (E %x% In + In %x% E)
}

# derivative of sigma0 with respect to sigma_zeta:
d_sigma0_sigma_zeta_var1 <- function(L,BetaStar,D2,...){
  L %*% BetaStar %*% D2
}

# Derivative of sigma_zeta to cholesky:
d_sigma_zeta_cholesky_var1 <- function(lowertri_zeta,L,C,In,...){
  d_sigma_cholesky(lowertri=lowertri_zeta,L=L,C=C,In=In)
}

# Derivative of sigma_zeta to precision:
d_sigma_zeta_kappa_var1 <- function(L,D2, sigma_zeta,...){
  d_sigma_kappa(L = L, D = D2, sigma = sigma_zeta)
}

# Derivative of sigma_zeta to ggm:
d_sigma_zeta_ggm_var1 <- function(L,delta_IminOinv_zeta,A,delta_zeta,Dstar,In,...){
  cbind(
    d_sigma_omega(L = L, delta_IminOinv = delta_IminOinv_zeta, A = A, delta = delta_zeta, Dstar = Dstar),
    d_sigma_delta(L = L,  delta_IminOinv = delta_IminOinv_zeta,In=In,delta=delta_zeta,A=A)
  )
}

## LAg 1
# Derivative of sigma1 with respect to beta:
d_sigma1_beta_var1 <- function(IkronBeta,D2,Jb,sigma,beta,In,...){
  n <- nrow(beta)
  sigma0 <- sigma[n + (1:n), n + (1:n)]
  as( IkronBeta %*% D2 %*% Jb + (sigma0 %x% In), "Matrix")
}


# Derivative of sigma1 with respect to omega:
d_sigma1_sigma_zeta_var1 <- function(IkronBeta,D2,Js,...){
  as(IkronBeta %*% D2 %*% Js, "Matrix")
}


# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_var1_group <- function(beta,P,zeta,...){
  dots <- list(...)
  
  # Number of variables:
  nvar <- nrow(beta) * 2
  
  # Number of nodes:
  nNode <- nvar / 2
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nvar + # Means
    nNode*(nNode+1)/2 + # Exogenous variances
    nNode^2 + # Beta
    nNode * (nNode+1) / 2 # Contemporaneous network and var-cov
  
  # Empty Jacobian:
  Jac <- matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaStarInds <- nvar + seq_len(nNode*(nNode+1)/2)
  sigma0Inds <- max(sigmaStarInds) + seq_len(nNode*(nNode+1)/2)
  sigma1Inds <- max(sigma0Inds) + seq_len(nNode^2)
  
  # Indices model:
  interceptInds <- 1:nvar
  exovarInds <-  nvar + seq_len(nNode*(nNode+1)/2)
  betaInds <- max(exovarInds) + seq_len(nNode^2)
  sigmazetaInds <-  max(betaInds) + seq_len(nNode*(nNode+1)/2)
  

  # fill intercept part:
  # Jac[meanInds,interceptInds] <- bdiag(d_mu_mu_var1(beta=beta,...),d_mu_mu_var1(beta=beta,...))
  Jac[meanInds,interceptInds] <- as.matrix(bdiag(d_mu_mu_var1_cpp(beta),d_mu_mu_var1_cpp(beta)))
  
  # Fill the exo var part:
  # Jac[sigmaStarInds,exovarInds] <- d_sigmastar_exo_cholesky_var1(...)
  Jac[sigmaStarInds,exovarInds] <- d_sigmastar_exo_cholesky_var1_cpp(In = dots$In, L = dots$L, C = dots$C, exo_cholesky = dots$exo_cholesky)
  
  # Fill sigma0 to beta part:
  # Jac[sigma0Inds,betaInds] <- Jb <- d_sigma0_beta_var1(...)
  Jac[sigma0Inds,betaInds] <- Jb <- d_sigma0_beta_var1_cpp(BetaStar = dots$BetaStar, In = dots$In, sigma = dots$sigma, C = dots$C, L = dots$L)
  
  # Fill sigma0 to sigma_zeta part:
  # Jac[sigma0Inds,sigmazetaInds] <- d_sigma0_sigma_zeta_var1(...)
  Jac[sigma0Inds,sigmazetaInds] <- d_sigma0_sigma_zeta_var1_cpp(L = dots$L, BetaStar = dots$BetaStar, D2 = dots$D2)
  # d_sigma0_sigma_zeta_var1_cpp
  
  # Augment:
  # if (zeta == "chol"){
  #   Jac[sigma0Inds,sigmazetaInds] <- as.matrix( Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_cholesky_var1(...) )
  # } else if (zeta == "prec"){
  #   Jac[sigma0Inds,sigmazetaInds] <- as.matrix( Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_kappa_var1(...))
  # } else if (zeta == "ggm"){
  #   Jac[sigma0Inds,sigmazetaInds] <- as.matrix(Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_ggm_var1(...))
  # }
  if (zeta == "chol"){
    Jac[sigma0Inds,sigmazetaInds] <-  Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_cholesky_var1_cpp(lowertri_zeta = dots$lowertri_zeta,L = dots$L,C = dots$C,In = dots$In)
  } else if (zeta == "prec"){
    Jac[sigma0Inds,sigmazetaInds] <-  Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_kappa_var1_cpp(L = dots$L,D2 = dots$D2,sigma_zeta = dots$sigma_zeta)
    
      # Jac[sigma0Inds,sigmazetaInds] <- as.matrix( Jac[sigma0Inds,sigmazetaInds] %*% as.matrix(d_sigma_zeta_kappa_var1(...)))
    
  } else if (zeta == "ggm"){
    Jac[sigma0Inds,sigmazetaInds] <- Jac[sigma0Inds,sigmazetaInds] %*% d_sigma_zeta_ggm_var1_cpp(L = dots$L, delta_IminOinv_zeta = dots$delta_IminOinv_zeta,A = dots$A, delta_zeta = dots$delta_zeta,Dstar = dots$Dstar,In = dots$In)
  }
  
  # Store:
  Js <- Jac[sigma0Inds,sigmazetaInds]
  
##
  # Fill sigma1 to beta part:
  # Jac[sigma1Inds,betaInds] <- as.matrix(d_sigma1_beta_var1(beta=beta,Jb=Jb,...))
  # 
  # # Fill sigma1 to sigma_zeta part:
  # Jac[sigma1Inds,sigmazetaInds] <- as.matrix(d_sigma1_sigma_zeta_var1(Js=Js,...))
  
  Jac[sigma1Inds,betaInds] <- d_sigma1_beta_var1_cpp(beta=beta,Jb=Jb,IkronBeta = dots$IkronBeta,D2 = dots$D2,sigma = dots$sigma,In = dots$In)

  # Fill sigma1 to sigma_zeta part:
  Jac[sigma1Inds,sigmazetaInds] <- d_sigma1_sigma_zeta_var1_cpp(Js=Js,IkronBeta = dots$IkronBeta, D2 = dots$D2)
  # 
  # Augment:
  # if (zeta == "chol"){
  #   Jac[sigma1Inds,sigmazetaInds] <- Jac[sigma1Inds,sigmazetaInds]  %*% d_sigma_zeta_cholesky_var1(...)
  # } else if (zeta == "prec"){
  #   Jac[sigma1Inds,sigmazetaInds] <- Jac[sigma1Inds,sigmazetaInds]  %*% d_sigma_zeta_kappa_var1(...)
  # } else if (zeta == "ggm"){
  #   Jac[sigma1Inds,sigmazetaInds] <- Jac[sigma1Inds,sigmazetaInds]  %*% d_sigma_zeta_ggm_var1(...)
  # }
  
  # Permute the matrix:
  Jac <- P %*% Jac

  # Make sparse if needed:
  Jac <- sparseordense(Jac)
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_var1 <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_var1_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




