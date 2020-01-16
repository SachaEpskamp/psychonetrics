# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_Ising_group <- function(omega,...){
  
  # Number of variables:
  nvar <- nrow(omega)
  
  # Jacobian is simply an identity matrix of the correct dimension:
  Jac <- Diagonal(nvar + (nvar*(nvar-1)/2) + 1)
    
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


