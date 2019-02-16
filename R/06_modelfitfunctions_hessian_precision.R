# Fit function per group:
hessian_precision_group <- function(S,kappa,means,mu,D,sigma,...){

  # Mean-mean part:
  hes_mean_mean <- 2*kappa
  
  # mean - net part:
  hes_mean_kappa <- -2 * (kronecker(t(means - mu),diag(length(means))) %*% D)
  
  # net - mean part:
  hes_kappa_mean <- -t(D) %*% (kronecker(diag(length(means)),means - mu) + kronecker(means - mu,diag(length(means))))
  
  # Check if equal:
  if (!all(hes_kappa_mean == t(hes_mean_kappa))){
    warning("Hessian is not symmetrical... Symmetrizing...")
    avg <- 0.5 * (hes_kappa_mean + t(hes_mean_kappa))
    hes_kappa_mean <- avg
    hes_mean_kappa <- t(avg)
  }
  
  # Network-network part:
  hes_kappa_kappa <- t(D) %*% kronecker(sigma, sigma) %*% D
  
  # Combine and return:
  cbind(
    rBind(hes_mean_mean,hes_kappa_mean),
    rBind(hes_mean_kappa,hes_kappa_kappa)
  )
}

# Fit function for the precision: -2n* log likelihood
hessian_precision <- function(x, model){
  return(as.matrix(hessian_general_kappa(x,model)))
  # # Prepare
  # prep <- prepare_precision(x, model)
  # 
  # # Fit function per group:
  # hessian_per_group <- lapply(prep$groupModels,do.call,what=hessian_precision_group)
  # 
  # # Make block matrix:
  # full_hessian <- Reduce("bdiag",hessian_per_group)
  # 
  # # Add model:
  # hes <- t(model@fitfunctions$extramatrices$M) %*% full_hessian %*% model@fitfunctions$extramatrices$M
  # 
  # # Return:
  # return(as.matrix(hes))
}