# Cholesky derivative:
d_sigma_cholesky <- function(lowertri,L,C,In,...){

  res <- L %*% ((In %x% In) + C) %*% ((lowertri %x% In) %*% t(L))
  
  as.matrix(res)
}

# Derivative of scaling matrix:
d_sigma_delta <- function(L,delta_IminOinv,In,A,delta,...){
  res <- L %*% (
    (delta_IminOinv%x% In) + 
      (In %x% delta_IminOinv)
  ) %*% A
  
  as.matrix(res)
}

# Derivative of network matrix:
d_sigma_omega <- function(L,delta_IminOinv,A,delta,Dstar,...){
  # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
  
  # delta_IminOinv <- delta %*% IminOinv
  res <- L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
  
  # all(a == b)
  as.matrix(res)
}

# Derivative of precision matrix:
d_sigma_kappa <- function(L,D,sigma,...){
  res <- - (L %*% (sigma %x% sigma) %*% D)
  as.matrix(res)
}

# Derivative of rho:
d_sigma_rho <- function(L,SD,A,delta,Dstar,...){
  res <- L %*% (SD %x% SD) %*% Dstar
  as.matrix(res)
}

# Derivative of SDs:
d_sigma_SD <- function(L,SD_IplusRho,In,A,...){
  res <- L %*% (
    (SD_IplusRho%x% In) + 
      (In %x% SD_IplusRho)
  ) %*% A
  as.matrix(res)
}

# Derivative of omega to covariances in the corinput = TRUE setting:
d_sigma_omega_corinput <- function(L,delta_IminOinv,A,delta,Dstar,IminOinv,In,...){
  # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
  
  # delta_IminOinv <- delta %*% IminOinv
  res <- L %*% (
    (delta_IminOinv %x% delta_IminOinv) -
      1/2 *  ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A %*% 
      Diagonal(x = diag(IminOinv)^(-1.5)) %*% t(A) %*% (IminOinv %x% IminOinv)
  ) %*% Dstar
  as.matrix(res)
}



# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_varcov_group <- function(cpp, sigma,y,corinput,meanstructure,tau,mu,...){
  dots <- list(...)
  
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
  Jac <- matrix(0, nobs, npars)
  
  if (meanstructure || nThresh > 0){
    # Fill mean part with diagonal:
    Jac[meanPart,meanPart] <- as.matrix(Diagonal(nMean_Thresh))
  }

  # Now fill the sigma part:
  if (y == "cov"){
    # Regular covs:
    Jac[varPart,varPartPars] <- as.matrix(Diagonal(nvar*(nvar+1)/2))
  } else if (y == "chol"){
  
    if (cpp){
      # Cholesky decomposition:
      Jac[varPart,varPartPars] <- d_sigma_cholesky_cpp(lowertri = dots$lowertri, L = dots$L, C = dots$C, In = dots$In)
    } else {
      # Cholesky decomposition:
      Jac[varPart,varPartPars] <- d_sigma_cholesky(...)      
    }

  } else if (y == "ggm"){
    
    # Gaussian graphical model:
    netPart <- meanstructure*max(meanPart) + nThresh + seq_len(nvar*(nvar-1)/2)
    scalingPart <- max(netPart) + seq_len(nvar)
    
    if (corinput){
      
      
      if (cpp){
        Jac[varPart,netPart] <- d_sigma_omega_corinput_cpp(L = dots$L, 
                                                           delta_IminOinv = dots$delta_IminOinv, 
                                                           A = dots$A, 
                                                           delta = dots$delta,
                                                           Dstar = dots$Dstar,
                                                           IminOinv = dots$IminOinv,
                                                           In = dots$In)
      } else {
        Jac[varPart,netPart] <- d_sigma_omega_corinput(...)
      }

      
    } else {
     if (cpp){
       Jac[varPart,netPart] <- d_sigma_omega_cpp(L = dots$L,delta_IminOinv = dots$delta_IminOinv,A = dots$A,delta = dots$delta,Dstar = dots$Dstar)
       Jac[varPart,scalingPart] <- d_sigma_delta_cpp(L = dots$L, delta_IminOinv = dots$delta_IminOinv, In = dots$In,A = dots$A)
     } else {
       Jac[varPart,netPart] <- d_sigma_omega(...)
       Jac[varPart,scalingPart] <- d_sigma_delta(...)       
     }

      
    }
    
  } else  if (y == "prec"){
    if (cpp){
      Jac[varPart,varPartPars] <- d_sigma_kappa_cpp(sigma=sigma,L = dots$L, D = dots$D)
    } else {
      Jac[varPart,varPartPars] <- d_sigma_kappa(sigma=sigma,...)      
    }

  } else if (y == "cor"){
    # Corelation matrix:
    corPart <- meanstructure*max(meanPart) + nThresh +  seq_len(nvar*(nvar-1)/2)
    
    if  (cpp){
      Jac[varPart,corPart] <- d_sigma_rho_cpp(L = dots$L, SD = dots$SD, A = dots$A, Dstar = dots$Dstar)
    } else {
      Jac[varPart,corPart] <- d_sigma_rho(...)      
    }

    
    if (!corinput){
      sdPart <- max(corPart) + seq_len(nvar)  
      
      if (cpp){
        Jac[varPart,sdPart] <- d_sigma_SD_cpp(L = dots$L, SD_IplusRho = dots$SD_IplusRho,In = dots$In, A = dots$A)
      } else {
        Jac[varPart,sdPart] <- d_sigma_SD(...)        
      }

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
    } else if (any(colSums(is.na(tau)) == nrow(tau))) stop("Mix of continuous and ordinal variables is not yet supported.")
  }

  # Make sparse if needed:
  Jac <- sparseordense(Jac)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_varcov <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_varcov_group)
  
  # FIXME: Computationally it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  sparseordense(Reduce("bdiag",d_per_group))
}


