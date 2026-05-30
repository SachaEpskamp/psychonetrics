# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_Ising_group <- function(omega,beta_model,log_beta,...){

  # Number of variables:
  nvar <- nrow(omega)

  # Distribution parameters phi = (tau, omega[lower.tri], delta, beta). The map
  # from model parameters theta to phi is the identity (both the Ising and
  # BlumeCapel models map directly onto the Spin distribution); for the Ising
  # model the delta entries are fixed to zero in the parameter table, so those
  # identity columns are simply never used. beta remains the last parameter.
  Jac <- Diagonal(nvar + (nvar*(nvar-1)/2) + nvar + 1)

  # Fix log beta:
  if (beta_model == "log_beta"){
    Jac[nrow(Jac),ncol(Jac)] <- exp(log_beta)
  }
    
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_Ising <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_Ising_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}


