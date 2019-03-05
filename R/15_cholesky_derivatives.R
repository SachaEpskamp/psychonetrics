
# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_cholesky_group <- function(mu,lowertri,L,C,In,...){
  # Number of variables:
  nvar <- length(mu)

  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # Indices:
  meanInds <- 1:nvar
  sigmaInds <- nvar + seq_len(nvar*(nvar+1)/2)
  
  # Empty Jacobian:
  Jac <- Matrix(0,nobs,nobs)

  # Diagonal for mu:
  Jac[meanInds,meanInds] <- Diagonal(n=nvar)
  
  # Cholesky derivative:
  Jac[sigmaInds,sigmaInds] <- L %*% (
    ( In %(x)% lowertri) %*% C %*% t(L) + 
      (lowertri %(x)% In) %*% t(L)
  )
  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_cholesky <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_cholesky_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}


