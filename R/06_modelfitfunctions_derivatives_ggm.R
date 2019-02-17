# Derivative of mu with respect to mu:
d_mu_mu_ggm <- function(mu,...){
  Diagonal(length(mu))
}

# Derivative of vech(kappa) with respect to vechs(omega):
d_kappa_omega_ggm <- function(L,delta,Dstar,...){
 -L %*% kronecker(delta, delta) %*% Dstar
}

# Derivative of vech(kappa) with respect to diag(delta):
d_kappa_delta_ggm <- function(L,delta,omega,E,...){
  I <- Diagonal(nrow(omega))
  L %*% (kronecker(delta %*% (I - omega), I) + kronecker(I, (I - omega)%*%delta))%*%E
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
    (nvar * (nvar-1))/2 + # edges
    nvar # scaling

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # fill mean part:
  Jac[1:nvar,1:nvar] <- d_mu_mu_ggm(...)
  
  # Fill network part:
  Jac[nvar + seq_len((nvar * (nvar+1))/2), nvar + seq_len((nvar * (nvar-1))/2)] <- d_kappa_omega_ggm(...)

  # Fill scaling part:
  indsx <- nvar  + seq_len((nvar * (nvar+1))/2)
  indsy <- nvar + ((nvar * (nvar-1))/2) + seq_len(nvar)
  Jac[indsx,indsy] <- d_kappa_delta_ggm(omega,...)
 
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
