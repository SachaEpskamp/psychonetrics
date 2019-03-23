
 # derivative of sigma0 (full vec) with respect to beta:
d_sigma0_beta_panelvar1 <- function(BetaStar,E,In,w,B_Sigma_mu,C,...){
  (t(w) %x% (In %x% In)) %*% 
    (t(BetaStar) %x% BetaStar) %*% 
    (E %x% In + In %x% E) - 
    BetaStar %*% (
      (In %x% B_Sigma_mu) %*% C + (B_Sigma_mu %x% In)
    )
}
# 
# # derivative of sigma0 with respect to sigma_zeta:
d_sigma0_sigma_zeta_panelvar1 <- function(BetaStar,D2,...){
  BetaStar %*% D2
}

## Lag k derivatives:
d_sigmak_beta_panelvar1 <- function(k,IkronBeta,Jb,allSigmas,Sigma_mu_kron_I,In,...){
  IkronBeta %*% Jb + (t(allSigmas[[k-1]]) %x% In) - Sigma_mu_kron_I
}

d_sigmak_sigma_zeta_panelvar1 <- function(IkronBeta,Jsigz,...){
  IkronBeta %*% Jsigz
}

d_sigmak_sigma_mu_panelvar1 <- function(D2,...){
  D2
}

### Augmentation parts:

# Derivative of sigma_zeta to cholesky:
d_sigma_zeta_cholesky_panelvar1 <- function(lowertri_zeta,L,C,In,...){
  d_sigma_cholesky(lowertri=lowertri_zeta,L=L,C=C,In=In)
}

# Derivative of sigma_zeta to precision:
d_sigma_zeta_kappa_panelvar1 <- function(L,D2, sigma_zeta,...){
  d_sigma_kappa(L = L, D = D2, sigma = sigma_zeta)
}

# Derivative of sigma_zeta to ggm:
d_sigma_zeta_ggm_panelvar1 <- function(L,IminOinv_zeta,A,delta_zeta,Dstar,In,...){
  cbind(
    d_sigma_omega(L = L, IminOinv = IminOinv_zeta, A = A, delta = delta_zeta, Dstar = Dstar),
    d_sigma_delta(L = L,  IminOinv = IminOinv_zeta,In=In,delta=delta_zeta,A=A)
  )
}


d_sigma_mu_cholesky_panelvar1 <- function(lowertri_mu,L,C,In,...){
  d_sigma_cholesky(lowertri=lowertri_mu,L=L,C=C,In=In)
}

# Derivative of sigma_mu to precision:
d_sigma_mu_kappa_panelvar1 <- function(L,D2, sigma_mu,...){
  d_sigma_kappa(L = L, D = D2, sigma = sigma_mu)
}

# Derivative of sigma_mu to ggm:
d_sigma_mu_ggm_panelvar1 <- function(L,IminOinv_mu,A,delta_mu,Dstar,In,...){
  cbind(
    d_sigma_omega(L = L, IminOinv = IminOinv_mu, A = A, delta = delta_mu, Dstar = Dstar),
    d_sigma_delta(L = L,  IminOinv = IminOinv_mu,In=In,delta=delta_mu,A=A)
  )
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_panelvar1_group <- function(...,design,P,contemporaneous,between,L){

  # Number of nodes:
  nNode <- nrow(design)
  
  # Number of times points:
  nTime <- ncol(design)
  
  # I need to construct the Jacobian for the following "observations:"
  nobs <- nNode + # Means
    (nNode * (nNode+1))/2 + # Variances
    nNode^2 * (nTime - 1) # lagged variances
  
  # total number of elements:
  nelement <- nNode + # Means
    nNode^2 + # Temporal effects
    nNode*(nNode+1)/2 + # Contemporaneous effects
    nNode*(nNode+1)/2 # Between-subject effects
  
  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nNode
  sigInds <- list(
    nNode + seq_len(nNode*(nNode+1)/2)
  )
  
  # For each lag:
  for (t in 2:nTime){
    sigInds[[t]] <- max(sigInds[[t-1]]) + seq_len(nNode^2)
  }
  
  # Indices model:
  interceptInds <- 1:nNode
  betaInds <- max(interceptInds) + seq_len(nNode^2)
  sigmazetaInds <-  max(betaInds) + seq_len(nNode*(nNode+1)/2)
  sigmamuInds <-  max(sigmazetaInds) + seq_len(nNode*(nNode+1)/2)
  
  # Augmentation parts:
  if (contemporaneous == "chol"){
    contAug <-  d_sigma_zeta_cholesky_panelvar1(L=L,...)
  } else if (contemporaneous == "prec"){
    contAug <-  d_sigma_zeta_kappa_panelvar1(L=L,...)
  } else if (contemporaneous == "ggm"){
    contAug <-  d_sigma_zeta_ggm_panelvar1(L=L,...)
  } else if (contemporaneous == "cov"){
    contAug <- Diagonal(nNode*(nNode+1)/2)
  }
  
  if (between == "chol"){
    betAug <-  d_sigma_mu_cholesky_panelvar1(L=L,...)
  } else if (between == "prec"){
    betAug <-  d_sigma_mu_kappa_panelvar1(L=L,...)
  } else if (between == "ggm"){
    betAug <-  d_sigma_mu_ggm_panelvar1(L=L,...)
  } else if (between == "cov"){
    betAug <- Diagonal(nNode*(nNode+1)/2)
  }

  # fill intercept part:
  Jac[meanInds,interceptInds] <- Diagonal(nNode)
  
  # Fill sigma0 to beta part:
  Jb <- d_sigma0_beta_panelvar1(...)
  Jac[sigInds[[1]],betaInds] <- L %*% Jb

  # Fill sigma0 to sigma_zeta part:
  Jsigz <- d_sigma0_sigma_zeta_panelvar1(...)
  Jac[sigInds[[1]],sigmazetaInds] <- L %*% Jsigz  %*% contAug
  
  # Fill sigma0 to sigma_mu part:
  Jac[sigInds[[1]],sigmamuInds] <- Diagonal(nNode*(nNode+1)/2) %*% betAug
  
  # For every further lag:
  for (t in 2:nTime){
    # Beta part:
    Jb <- d_sigmak_beta_panelvar1(k=t,Jb=Jb,...)
    Jac[sigInds[[t]],betaInds] <- Jb
    
    # sigma_zeta part:
    Jsigz <- d_sigmak_sigma_zeta_panelvar1(Jsigz=Jsigz,...)
    Jac[sigInds[[t]],sigmazetaInds] <- Jsigz %*% contAug
    
    # sigma_omega part:
    Jac[sigInds[[t]],sigmamuInds] <- d_sigmak_sigma_mu_panelvar1(...) %*% betAug
  }
  

  # Permute the matrix:
  Jac <- P %*% Jac
  

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_panelvar1 <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_panelvar1_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




