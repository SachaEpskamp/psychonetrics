# Derivative of tau with respect to mu:
d_mu_mu_gvar <- function(mu,...){
  Diagonal(length(mu))
}

# Derivative of exogenous variances part
d_sigmastar_sigmastar_gvar <- function(exogenous_sigma,...){
  n <- nrow(exogenous_sigma)
  Diagonal(n*(n+1)/2)
}

# derivative of sigma0 with respect to beta:
d_sigma0_beta_gvar <- function(L_betaStar,E,In,sigmaZetaVec,BetaStar,...){
  ((t(sigmaZetaVec) %*% t(BetaStar)) %(x)% L_betaStar) %*% (E %(x)% In + In %(x)% E)
}

# derivative of sigma0 with respect to omega:
d_sigma0_omega_gvar <- function(L,BetaStar,delta_zeta,OmegaStar,Dstar,...){
  L %*% BetaStar %*% (delta_zeta %(x)% delta_zeta) %*% (OmegaStar %(x)% OmegaStar) %*% Dstar
}

# derivative of sigma0 with respect to delta:
d_sigma0_delta_gvar <- function(L,BetaStar,DeltaOmegaStar,A,In,...){
  L %*% BetaStar %*% (DeltaOmegaStar %(x)% In + In %(x)% DeltaOmegaStar) %*% A
}

# Derivative of sigma1 with respect to beta:
d_sigma1_beta_gvar <- function(IkronBeta,D2,Jb,sigma,beta,In,...){
  n <- nrow(beta)
  sigma0 <- sigma[n + (1:n), n + (1:n)]
  IkronBeta %*% D2 %*% Jb + (sigma0 %(x)% In)
}

# Derivative of sigma1 with respect to omega:
d_sigma1_omega_gvar <- function(IkronBeta,D2,Jo,...){
  IkronBeta %*% D2 %*% Jo
}

# Derivative of sigma1 with respect to delta:
d_sigma1_delta_gvar <- function(IkronBeta,D2,Jd,...){
  IkronBeta %*% D2 %*% Jd
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_gvar_group <- function(beta,P,...){
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
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaStarInds <- nvar + seq_len(nNode*(nNode+1)/2)
  sigma0Inds <- max(sigmaStarInds) + seq_len(nNode*(nNode+1)/2)
  sigma1Inds <- max(sigma0Inds) + seq_len(nNode^2)
  
  # Indices model:
  interceptInds <- 1:nvar
  exovarInds <-  nvar + seq_len(nNode*(nNode+1)/2)
  betaInds <- max(exovarInds) + seq_len(nNode^2)
  omegaInds <-  max(betaInds) + seq_len(nNode*(nNode-1)/2)
  deltaInds <- max(omegaInds) + seq_len(nNode)

  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_mu_gvar(...)
  
  # Fill the exo var part:
  Jac[sigmaStarInds,exovarInds] <- d_sigmastar_sigmastar_gvar(...)
  
  # Fill sigma0 to beta part:
  Jac[sigma0Inds,betaInds] <- Jb <- d_sigma0_beta_gvar(...)
  
  # Fill sigma0 to omega part:
  Jac[sigma0Inds,omegaInds] <- Jo <- d_sigma0_omega_gvar(...)
  
  # Fill sigma0 to delta part:
  Jac[sigma0Inds,deltaInds] <- Jd <- d_sigma0_delta_gvar(...)
  
##
  # Fill sigma1 to beta part:
  Jac[sigma1Inds,betaInds] <- d_sigma1_beta_gvar(beta=beta,Jb=Jb,...)
  
  # Fill sigma1 to omega part:
  Jac[sigma1Inds,omegaInds] <- d_sigma1_omega_gvar(Jo=Jo,...)
  
  # Fill sigma1 to delta part:
  Jac[sigma1Inds,deltaInds] <- d_sigma1_delta_gvar(Jd=Jd,...)
  
  # Permute the matrix:
  Jac <- P %*% Jac

  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_gvar <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_gvar_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




