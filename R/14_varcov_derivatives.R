
# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_varcov_group <- function(sigma,...){
  # Number of variables:
  nvar <- nrow(sigma)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # Diagonal Jacobian:
  Jac <- Diagonal(n = nobs)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_varcov <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_varcov_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}


