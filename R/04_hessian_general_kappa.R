# Fit function for the precision: -2n* log likelihood
hessian_general_kappa <- function(x, model){
  # Prepare
  prep <- prepare_precision(x, model)

  # Fit function per group:
  hessian_per_group <- lapply(prep$groupModels,do.call,what=hessian_precision_group)
  
  # Make block matrix:
  full_hessian <- Reduce("bdiag",hessian_per_group)
  
  # Add model:
  hes <- t(model@fitfunctions$extramatrices$M) %*% full_hessian %*% model@fitfunctions$extramatrices$M
  
  # Return:
  return(hes)
}