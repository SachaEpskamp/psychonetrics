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

# Derivative of rho:
d_sigma_rho <- function(L,SD,A,delta,Dstar,...){
  L %*% (SD %x% SD) %*% Dstar
}

# Derivative of SDs:
d_sigma_SD <- function(L,SD_IplusRho,In,A,...){
  L %*% (
    (SD_IplusRho%x% In) + 
      (In %x% SD_IplusRho)
  ) %*% A
}

# Derivative of omega to covariances in the corinput = TRUE setting:
d_sigma_omega_corinput <- function(L,delta_IminOinv,A,delta,Dstar,IminOinv,In,...){
  # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
  
  # delta_IminOinv <- delta %*% IminOinv
  L %*% (
    (delta_IminOinv %x% delta_IminOinv) -
      1/2 *  ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A %*% 
      Diagonal(x = diag(IminOinv)^(-1.5)) %*% t(A) %*% (IminOinv %x% IminOinv)
  ) %*% Dstar
}



# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_varcov_group <- function(sigma,y,corinput,meanstructure,tau,mu,...){
  
  # Number of variables:
  nvar <- nrow(sigma)
  
  if (missing(tau)){
    tau <- matrix(NA,1,nvar)
  }
  

  # Number of means/thresholds:
  nMean_Thresh <- sum(!is.na(tau)) + sum(!is.na(mu))
  
  nThresh <- sum(!is.na(tau))
  
  # Number of observations:
  nobs <- nMean_Thresh + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # Number of parameters is less if corinput is used or if meanstructure is ignored:
  npars <- nobs - corinput * nvar - (!meanstructure) * sum(!is.na(mu))
  
  # Mean part:
  meanPart <- seq_len(nMean_Thresh)
  
  # Variance part:
  varPart <- max(meanPart) + seq_len(nvar*(nvar+1)/2)    
  
  # Var part for parameters:
  varPartPars <- meanstructure * max(meanPart) + nThresh +  seq_len(nvar*(nvar+1)/2)    
  
  # Empty Jacobian:
  Jac <- Matrix(0, nobs, npars, sparse = FALSE)
  
  if (meanstructure || nThresh > 0){
    # Fill mean part with diagonal:
    Jac[meanPart,meanPart] <- Diagonal(nMean_Thresh)    
  }
  
  
  # Now fill the sigma part:
  if (y == "cov"){
    # Regular covs:
    Jac[varPart,varPartPars] <- Diagonal(nvar*(nvar+1)/2)
  } else if (y == "chol"){
    # Cholesky decomposition:
    Jac[varPart,varPartPars] <- d_sigma_cholesky(...)
  } else if (y == "ggm"){
    
    # Gaussian graphical model:
    netPart <- meanstructure*max(meanPart) + nThresh + seq_len(nvar*(nvar-1)/2)
    scalingPart <- max(netPart) + seq_len(nvar)
    
    if (corinput){
      
      Jac[varPart,netPart] <- d_sigma_omega_corinput(...)
      
    } else {
     
      Jac[varPart,netPart] <- d_sigma_omega(...)
      Jac[varPart,scalingPart] <- d_sigma_delta(...)
      
    }
    
  } else  if (y == "prec"){
    Jac[varPart,varPartPars] <- d_sigma_kappa(sigma=sigma,...)
  } else if (y == "cor"){
    # Corelation matrix:
    corPart <- meanstructure*max(meanPart) + nThresh +  seq_len(nvar*(nvar-1)/2)
    Jac[varPart,corPart] <- d_sigma_rho(...)
    
    if (!corinput){
      sdPart <- max(corPart) + seq_len(nvar)  
      Jac[varPart,sdPart] <- d_sigma_SD(...)
    }
  }
  
  # Cut out the rows not needed
  # FIXME: Nicer to not have to compute these in the first place...
  if (corinput){
    keep <- c(rep(TRUE,nMean_Thresh),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Jac <- Jac[keep,]
  }
  if (!meanstructure){
    if (all(is.na(tau))){
      Jac <- Jac[-(seq_len(nvar)), ] 
    } else if (any(is.na(tau))) stop("Mix of continuous and ordinal variables is not yet supported.")
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


