# Derivative of mu with respect to mu:
d_mu_mu_ggm <- function(mu,...){
  Diagonal(length(mu))
}


# Derivative of scaling matrix:
d_sigma_delta_ggm <- function(L,IminOinv,In,A,delta,...){
  L %*% (
    ((delta %*% IminOinv) %(x)% In) + 
      (In %(x)% (delta %*% IminOinv))
  ) %*% A
}

# Derivative of network matrix:
d_sigma_omega_ggm <- function(L,IminOinv,A,delta,Dstar,...){
  L %*% (delta %(x)% delta) %*% (IminOinv %(x)% IminOinv) %*% Dstar
}


# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_ggm_group <- function(omega,...){
  # Number of variables:
  nvar <- nrow(omega)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nvar + # Means
    nvar*(nvar + 1)/2 # scaling and network

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaInds <- nvar + seq_len(nvar*(nvar+1)/2)
  
  # Indices model:
  interceptInds <- 1:nvar
  omegainds <- max(interceptInds) + seq_len(nvar*(nvar-1)/2)
  deltainds <- max(omegainds) + seq_len(nvar)


  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_mu_ggm(...)
  
    # Fill network part:
  Jac[sigmaInds,omegainds] <- d_sigma_omega_ggm(lambda=lambda,...)
  
  # Fill scaling part:
  Jac[sigmaInds,deltainds] <- d_sigma_delta_ggm(lambda=lambda,...)
  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_ggm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_ggm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}


