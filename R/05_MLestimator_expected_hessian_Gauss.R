expected_hessian_Gaussian_group_meanPart <- function(kappa,...){
  2*kappa
}

expected_hessian_Gaussian_group_varPart <- function(kappa,D,Drawts,...){
  # if (nrow(Drawts) != ncol(Drawts)){
  #   D <- Diagonal(n = nrow(kappa)^2)
  # } 
  # Make kappa sparse:
  kappa[abs(kappa) < sqrt(.Machine$double.eps)] <- 0
  kappa <- as(kappa, "Matrix")
  
  t(D) %*% (kappa %x% kappa) %*% D
}

# Per group:
expected_hessian_Gaussian_group <- function(...,mu,sigma,corinput = FALSE,meanstructure = TRUE){
  # if (missing(Drawts)){
  #   Drawts <- Diagonal(NROW(mu) +  nrow(sigma) * ( nrow(sigma)+1) / 2)
  # }
  
  # Variance part:
  varPart <- expected_hessian_Gaussian_group_varPart(...)
  
  # Cut out variances if needed:
  if (corinput){
    keep <- diag(ncol(sigma))[lower.tri(diag(ncol(sigma)),diag=TRUE)] != 1
    varPart <- as(varPart[keep,keep, drop=FALSE], "Matrix")
  }
  
  # Cut out means if needed:
  if (meanstructure){
    # Mean part:
    meanPart <- expected_hessian_Gaussian_group_meanPart(...)
    
    # Put int output:
    Out <- bdiag(meanPart, varPart) 
  } else {
    Out <- varPart
  }
  
  

  # Return as block matrix:
  return(Out)
}

# Total:
expected_hessian_Gaussian <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_Gaussian_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",exph_per_group)
}
