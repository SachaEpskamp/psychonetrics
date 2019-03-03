# # Derivative of tau with respect to mu:
# d_mu_mu_gvar <- function(beta,...){
#   Diagonal(nrow(beta))
# }
# 
# # Derivative of exogenous variances part
# d_sigmastar_sigmastar_gvar <- function(exogenous_sigma,...){
#   n <- nrow(exogenous_sigma)
#   Diagonal(n*(n+1)/2)
# }
# 
# # derivative of sigma0 with respect to beta:
# d_sigma0_beta_gvar <- function(L_betaStar,E,In,sigmaZetaVec,BetaStar,...){
#   ((t(sigmaZetaVec) %*% t(BetaStar)) %(x)% L_betaStar) %*% (E %(x)% In + In %(x)% E)
# }
# 
# # derivative of sigma0 with respect to omega:
# d_sigma0_omega_gvar <- function(L,BetaStar,delta_zeta,OmegaStar,Dstar,...){
#   L %*% BetaStar %*% (delta_zeta %(x)% delta_zeta) %*% (OmegaStar %(x)% OmegaStar) %*% Dstar
# }
# 
# # derivative of sigma0 with respect to delta:
# d_sigma0_delta_gvar <- function(L,BetaStar,DeltaOmegaStar,A,In,...){
#   L %*% BetaStar %*% (DeltaOmegaStar %(x)% In + In %(x)% DeltaOmegaStar) %*% A
# }
# 
# # Derivative of sigma1 with respect to beta:
# d_sigma1_beta_gvar <- function(IkronBeta,D2,Jb,sigma,beta,In,...){
#   n <- nrow(beta)
#   sigma0 <- sigma[n + (1:n), n + (1:n)]
#   IkronBeta %*% D2 %*% Jb + (sigma0 %(x)% In)
# }
# 
# # Derivative of sigma1 with respect to omega:
# d_sigma1_omega_gvar <- function(IkronBeta,D2,Jo,...){
#   IkronBeta %*% D2 %*% Jo
# }
# 
# # Derivative of sigma1 with respect to delta:
# d_sigma1_delta_gvar <- function(IkronBeta,D2,Jd,...){
#   IkronBeta %*% D2 %*% Jd
# }

# Derivative of sigmak with respect to beta:
d_sigmak_beta_gvar <- function(IkronBeta,D2,Jb,sigma,beta,In,k,...){
  n <- nrow(beta)
  sigma0 <- sigma[n + (1:n), n + (1:n)]
  (In %(x)% (beta %^% k)) %*% D2 %*% Jb + Reduce("+",lapply(1:k,function(i){
    ((sigma0 %*% (t(beta)%^%(i-1))) %(x)% (beta%^%(k-i)))
  }))
}

# Derivative of sigmak with respect to omega:
d_sigmak_omega_gvar <- function(IkronBeta,D2,Jo,k,In,beta,...){
  (In %(x)% (beta %^% k)) %*%  D2 %*% Jo
}

# Derivative of sigmak with respect to delta:
d_sigmak_delta_gvar <- function(IkronBeta,D2,Jd,k,In,beta,...){
  (In %(x)% (beta %^% k)) %*%  D2 %*% Jd
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_gvar_group_rawts <- function(beta,P,mis,...){
  

  # Number of nodes:
  nNode <- ncol(mis)
  
  # Number of time points:
  nTime <- nrow(mis)
  
  # lags:
  nLag <- nTime-1
  
  # Number of distributional parameters:
  nobs <- 
    nNode + # Means
    nNode * (nNode+1)/2 + # stationary variances
    nLag * nNode^2 # lag-k varianves
    
  # total number of elements:
  nelement <- nNode + # Means
    nNode^2 + # Beta
    nNode * (nNode+1) / 2 # Contemporaneous network and var-cov
  
  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nNode
  sigma0Inds <- max(meanInds) + seq_len(nNode*(nNode+1)/2)
  curMax <- max(sigma0Inds)
  sigmaInds <- list()
  for (i in 1:nLag){
    sigmaInds[[i]] <- curMax+ seq_len(nNode^2)
    curMax <- max(sigmaInds[[i]])
  }

  # Indices model:
  interceptInds <- 1:nNode
  betaInds <- max(interceptInds) + seq_len(nNode^2)
  omegaInds <-  max(betaInds) + seq_len(nNode*(nNode-1)/2)
  deltaInds <- max(omegaInds) + seq_len(nNode)

  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_mu_gvar(beta=beta,...)

  # Fill sigma0 to beta part:
  Jac[sigma0Inds,betaInds] <- Jb <- d_sigma0_beta_gvar(...)
  
  # Fill sigma0 to omega part:
  Jac[sigma0Inds,omegaInds] <- Jo <- d_sigma0_omega_gvar(...)
  
  # Fill sigma0 to delta part:
  Jac[sigma0Inds,deltaInds] <- Jd <- d_sigma0_delta_gvar(...)
  
##
  # for (i in 1:nLag){
    for (i in 1:4){ # FIXME: Let's try with lag 4 max
    # Fill sigma1 to beta part:
    Jac[sigmaInds[[i]],betaInds] <- d_sigmak_beta_gvar(beta=beta,Jb=Jb,k=i,...)
    
    # Fill sigma1 to omega part:
    Jac[sigmaInds[[i]],omegaInds] <- d_sigmak_omega_gvar(beta=beta,Jo=Jo,k=i,...)
    
    # Fill sigma1 to delta part:
    Jac[sigmaInds[[i]],deltaInds] <- d_sigmak_delta_gvar(beta=beta,Jd=Jd,k=i,...)
  }
  
  
 
  
  # Permute the matrix:
  # Jac <- P %*% Jac

  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_gvar_rawts <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_gvar_group_rawts)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




