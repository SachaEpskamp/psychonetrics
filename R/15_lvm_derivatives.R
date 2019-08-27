# Derivative of nu with respect to mu:
d_mu_nu_lvm <- function(nu,...){
  Diagonal(length(nu))
}

# derivative of latent intecepts:
d_mu_nu_eta_lvm <- function(Lambda_BetaStar,...){
  Lambda_BetaStar
}

# Derivative of factor loadings to means:
d_mu_lambda_lvm <- function(nu_eta,BetaStar,In,...){
  (t(nu_eta) %*% t(BetaStar)) %x% In
}

# Derivative of beta to means:
d_mu_beta_lvm <- function(nu_eta,lambda,tBetakronBeta,...){
  (t(nu_eta) %x% lambda) %*% tBetakronBeta
}


# Derivative of factor loadings to vars:
d_sigma_lambda_lvm <- function(L,Lambda_BetaStar,Betasta_sigmaZeta,In,C,...){
  L %*% (
    ((Lambda_BetaStar %*% t(Betasta_sigmaZeta)) %x% In) + 
      (In %x% (Lambda_BetaStar %*% t(Betasta_sigmaZeta)))%*%C
  )
}

# Derivative of beta to vars:
d_sigma_beta_lvm <- function(L, lambda, Betasta_sigmaZeta, Cbeta,Inlatent,tBetakronBeta, ... ){
  L %*% (lambda %x% lambda) %*% (
    (Betasta_sigmaZeta %x% Inlatent) + 
      (Inlatent %x% Betasta_sigmaZeta) %*% Cbeta
  ) %*% tBetakronBeta
}

# Derivative of latent variance-covariance matrix:
d_sigma_sigma_zeta_lvm <- function(L,Lambda_BetaStar,Deta,...){
  L %*% (Lambda_BetaStar %x% Lambda_BetaStar) %*% Deta
}

# Derivative of sigma_zeta to cholesky:
d_sigma_zeta_cholesky_lvm <- function(lowertri_zeta,L_eta,Cbeta,Inlatent,...){
  d_sigma_cholesky(lowertri=lowertri_zeta,L=L_eta,C=Cbeta,In=Inlatent)
}

# Derivative of sigma_zeta to precision:
d_sigma_zeta_kappa_lvm <- function(L_eta,Deta,sigma_zeta,...){
  d_sigma_kappa(L = L_eta, D = Deta, sigma = sigma_zeta)
}

# Derivative of sigma_zeta to ggm:
d_sigma_zeta_ggm_lvm <- function(L_eta,delta_IminOinv_zeta,Aeta,delta_zeta,Dstar_eta,Inlatent,...){
  cbind(
    d_sigma_omega(L = L_eta, delta_IminOinv = delta_IminOinv_zeta, A = Aeta, delta = delta_zeta, Dstar = Dstar_eta),
    d_sigma_delta(L = L_eta,  delta_IminOinv = delta_IminOinv_zeta,In=Inlatent,delta=delta_zeta,A=Aeta)
  )
}


# Residual vars:

# Derivative of sigma_epsilon to cholesky:
d_sigma_epsilon_cholesky_lvm <- function(lowertri_epsilon,L,C_chol,In,...){
  d_sigma_cholesky(lowertri=lowertri_epsilon,L=L,C=C_chol,In=In)
}

# Derivative of sigma_epsilon to precision:
d_sigma_epsilon_kappa_lvm <- function(L,D,sigma_epsilon,...){
  d_sigma_kappa(L = L, D = D, sigma = sigma_epsilon)
}

# Derivative of sigma_epsilon to ggm:
d_sigma_epsilon_ggm_lvm <- function(L,delta_IminOinv_epsilon,A,delta_epsilon,Dstar,In,...){
  cbind(
    d_sigma_omega(L = L, delta_IminOinv = delta_IminOinv_epsilon, A = A, delta = delta_epsilon, Dstar = Dstar),
    d_sigma_delta(L = L,  delta_IminOinv = delta_IminOinv_epsilon,In=In,delta=delta_epsilon,A=A)
  )
}



# 
# 
# # Derivative of residual network matrix:
# d_sigma_omega_epsilon_lvm <- function(L,delta_epsilon,OmegaStar,Dstar,...){
#   L %*% (delta_epsilon %x% delta_epsilon) %*% (OmegaStar %x% OmegaStar) %*% Dstar
# }
# 
# # Derivative residual scalign:
# d_sigma_delta_epsilon_lvm <- function(L,delta_epsilon,OmegaStar,In,A,...){
#   DO <- delta_epsilon %*% OmegaStar
#   L %*% (
#     (DO %x% In) + (In %x% DO)
#   ) %*% A
# }

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_lvm_group <- function(lambda,latent,residual,...){
  # Number of variables:
  nvar <- nrow(lambda)
  
  # Number of latents:
  nlat <- ncol(lambda)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nvar + # Means
    nlat + # Latent intercpets
    nvar * nlat + # factor loadings
    nlat^2 + # beta elements
    nlat*(nlat + 1)/2 + # Latent variances and covariances
    nvar *( nvar + 1)/2 # Residual network and scaling

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement, sparse = FALSE)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaInds <- nvar + seq_len(nvar*(nvar+1)/2)
  
  # Indices model:
  interceptInds <- 1:nvar
  nuetaInds <- nvar + seq_len(nlat)
  lambdaInds <- max(nuetaInds) + seq_len(nlat*nvar)
  betaInds <- max(lambdaInds) + seq_len(nlat^2)
  sigmazetaInds <- max(betaInds) + seq_len(nlat*(nlat+1)/2)
  sigmaepsilonInds <- max(sigmazetaInds) + seq_len(nvar*(nvar+1)/2)

  
  
  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_nu_lvm(...)
  
  # Fill latent intercept part:
  Jac[meanInds,nuetaInds] <- d_mu_nu_eta_lvm(...)
  
  # Fill factor loading parts:
  Jac[meanInds,lambdaInds] <- d_mu_lambda_lvm(...)
  Jac[sigmaInds,lambdaInds] <- d_sigma_lambda_lvm(...)
  
  # Fill the beta parts:
  Jac[meanInds,betaInds] <- d_mu_beta_lvm(lambda=lambda,...)
  Jac[sigmaInds,betaInds] <- d_sigma_beta_lvm(lambda=lambda,...)

  # Fill latent variances part:
  Jac[sigmaInds,sigmazetaInds] <- d_sigma_sigma_zeta_lvm(...)
  
  if (latent == "chol"){
    Jac[sigmaInds,sigmazetaInds] <-Jac[sigmaInds,sigmazetaInds] %*% d_sigma_zeta_cholesky_lvm(...)
  } else if (latent == "prec"){
    Jac[sigmaInds,sigmazetaInds] <-Jac[sigmaInds,sigmazetaInds] %*% d_sigma_zeta_kappa_lvm(...)
  } else if (latent == "ggm"){
    Jac[sigmaInds,sigmazetaInds] <- Jac[sigmaInds,sigmazetaInds] %*% d_sigma_zeta_ggm_lvm(...)
  }
  
  # Residual variances:
  if (residual == "cov"){
    Jac[sigmaInds,sigmaepsilonInds] <- Diagonal(nvar*(nvar+1)/2)  
  } else if (residual == "chol"){
    Jac[sigmaInds,sigmaepsilonInds] <- d_sigma_epsilon_cholesky_lvm(...)
  } else if (residual == "prec"){
    Jac[sigmaInds,sigmaepsilonInds] <-  d_sigma_epsilon_kappa_lvm(...)
  } else if (residual == "ggm"){
    Jac[sigmaInds,sigmaepsilonInds] <-  d_sigma_epsilon_ggm_lvm(...)
  }
  
  
  # # Fill residual network part:
  # Jac[sigmaInds,omegainds] <- d_sigma_omega_epsilon_lvm(...)
  # 
  # # Fill residual scaling part:
  # Jac[sigmaInds,deltainds] <- d_sigma_delta_epsilon_lvm(...)
  # Make sparse if needed:
  Jac <- as(Jac, "Matrix")
  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_lvm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_lvm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




