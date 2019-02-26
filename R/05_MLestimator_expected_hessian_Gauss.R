expected_hessian_Gaussian_group_meanPart <- function(kappa,...){
  2*kappa
}

expected_hessian_Gaussian_group_varPart <- function(kappa,D,...){
  t(D) %*% (kappa %(x)% kappa) %*% D
}

# Per group:
expected_hessian_Gaussian_group <- function(...){
  
  # Mean part:
  meanPart <- expected_hessian_Gaussian_group_meanPart(...)
 
  # Variance part:
  varPart <- expected_hessian_Gaussian_group_varPart(...)
  
  # Return as block matrix:
  bdiag(meanPart, varPart)
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
