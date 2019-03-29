# Cholesky derivative:
d_sigma_cholesky <- function(lowertri,L,C,In,...){

  # library(microbenchmark)
  # microbenchmark(
  #   L %*% (
  #     ( In %x% lowertri) %*% C %*% t(L) + 
  #       (lowertri %x% In) %*% t(L)
  #   ),
  #   L %*% ((In %x% In) + C) %*% ((lowertri %x% In) %*% t(L))
  # )
 # a <- L %*% (
 #    ( In %x% lowertri) %*% C %*% t(L) + 
 #      (lowertri %x% In) %*% t(L)
 #  )
 
L %*% ((In %x% In) + C) %*% ((lowertri %x% In) %*% t(L))
}

# Derivative of scaling matrix:
d_sigma_delta <- function(L,delta_IminOinv,In,A,delta,...){
  L %*% (
    (delta_IminOinv%x% In) + 
      (In %x% delta_IminOinv)
  ) %*% A
}

# Derivative of network matrix:
d_sigma_omega <- function(L,delta_IminOinv,A,delta,Dstar,...){
  # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
  
  # delta_IminOinv <- delta %*% IminOinv
  L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
  
  # all(a == b)
}

# Derivative of precision matrix:
d_sigma_kappa <- function(L,D,sigma,...){
  - L %*% (sigma %x% sigma) %*% D
}


# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_varcov_group <- function(sigma,y,...){
  # Number of variables:
  nvar <- nrow(sigma)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # Mean part:
  meanPart <- seq_len(nvar)
  
  # Variance part:
  varPart <- max(meanPart) + seq_len(nvar*(nvar+1)/2)
  
  # Empty Jacobian:
  Jac <- Matrix(0, nobs, nobs, sparse = FALSE)

  # Fill mean part with diagonal:
  Jac[meanPart,meanPart] <- Diagonal(nvar)
  
  # Now fill the sigma part:
  if (y == "cov"){
    # Regular covs:
    Jac[varPart,varPart] <- Diagonal(nvar*(nvar+1)/2)
  } else if (y == "chol"){
    # Cholesky decomposition:
    Jac[varPart,varPart] <- d_sigma_cholesky(...)
  } else if (y == "ggm"){
    # Gaussian graphical model:
    netPart <- max(meanPart) + seq_len(nvar*(nvar-1)/2)
    scalingPart <- max(netPart) + seq_len(nvar)
    
    
    Jac[varPart,netPart] <- d_sigma_omega(...)
    Jac[varPart,scalingPart] <- d_sigma_delta(...)
  } else  if (y == "prec"){
    Jac[varPart,varPart] <- d_sigma_kappa(sigma=sigma,...)
  }
  
  # Make sparse if needed:
  Jac <- as(Jac, "Matrix")
 
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


