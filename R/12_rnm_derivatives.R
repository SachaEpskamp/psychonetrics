# Derivative of tau with respect to mu:
d_mu_tau_rnm <- function(tau,...){
  Diagonal(length(tau))
}

# derivative of latent intecepts:
d_mu_tau_eta_rnm <- function(Lambda_BetaStar,...){
  Lambda_BetaStar
}

# Derivative of factor loadings to means:
d_mu_lambda_rnm <- function(tau_eta,BetaStar,In,...){
  (t(tau_eta) %*% t(BetaStar)) %(x)% In
}

# Derivative of beta to means:
d_mu_beta_rnm <- function(tau_eta,lambda,tBetakronBeta,...){
  (t(tau_eta) %(x)% lambda) %*% tBetakronBeta
}


# Derivative of factor loadings to vars:
d_sigma_lambda_rnm <- function(L,Lambda_BetaStar,Betasta_sigmaZeta,In,C,...){
  L %*% (
    ((Lambda_BetaStar %*% t(Betasta_sigmaZeta)) %(x)% In) + 
      (In %(x)% (Lambda_BetaStar %*% t(Betasta_sigmaZeta)))%*%C
  )
}

# Derivative of beta to vars:
d_sigma_beta_rnm <- function(L, lambda, Betasta_sigmaZeta, Cbeta,Inlatent,tBetakronBeta, ... ){
  L %*% (lambda %(x)% lambda) %*% (
    (Betasta_sigmaZeta %(x)% Inlatent) + 
      (Inlatent %(x)% Betasta_sigmaZeta) %*% Cbeta
  ) %*% tBetakronBeta
}

# Derivative of latent variance-covariance matrix:
d_sigma_sigma_zeta_rnm <- function(L,Lambda_BetaStar,Deta,...){
  L %*% (Lambda_BetaStar %(x)% Lambda_BetaStar) %*% Deta
}

# Derivative of residual network matrix:
d_sigma_omega_epsilon_rnm <- function(L,delta_epsilon,OmegaStar,Dstar,...){
  L %*% (delta_epsilon %(x)% delta_epsilon) %*% (OmegaStar %(x)% OmegaStar) %*% Dstar
}

# Derivative residual scalign:
d_sigma_delta_epsilon_rnm <- function(L,delta_epsilon,OmegaStar,In,A,...){
  DO <- delta_epsilon %*% OmegaStar
  L %*% (
    (DO %(x)% In) + (In %(x)% DO)
  ) %*% A
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_rnm_group <- function(lambda,...){
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
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaInds <- nvar + seq_len(nvar*(nvar+1)/2)
  
  # Indices model:
  interceptInds <- 1:nvar
  tauetaInds <- nvar + seq_len(nlat)
  lambdaInds <- max(tauetaInds) + seq_len(nlat*nvar)
  betaInds <- max(lambdaInds) + seq_len(nlat^2)
  sigmazetaInds <- max(betaInds) + seq_len(nlat*(nlat+1)/2)
  omegainds <- max(sigmazetaInds) + seq_len(nvar*(nvar-1)/2)
  deltainds <- max(omegainds) + seq_len(nvar)

  
  
  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_tau_rnm(...)
  
  # Fill latent intercept part:
  Jac[meanInds,tauetaInds] <- d_mu_tau_eta_rnm(...)
  
  # Fill factor loading parts:
  Jac[meanInds,lambdaInds] <- d_mu_lambda_rnm(...)
  Jac[sigmaInds,lambdaInds] <- d_sigma_lambda_rnm(...)
  
  # Fill the beta parts:
  Jac[meanInds,betaInds] <- d_mu_beta_rnm(lambda=lambda,...)
  Jac[sigmaInds,betaInds] <- d_sigma_beta_rnm(lambda=lambda,...)
  
  # Fill latent variances part:
  Jac[sigmaInds,sigmazetaInds] <- d_sigma_sigma_zeta_rnm(...)
  
  # Fill residual network part:
  Jac[sigmaInds,omegainds] <- d_sigma_omega_epsilon_rnm(...)
  
  # Fill residual scaling part:
  Jac[sigmaInds,deltainds] <- d_sigma_delta_epsilon_rnm(...)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_rnm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_rnm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}




