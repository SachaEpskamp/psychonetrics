# jacobian function per group:
hessian_gaussian_kappa_group <- function(S,kappa,means,mu,D,sigma,...){
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

# Now for all groups:
hessian_gaussian_kappa <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  h_per_group <- lapply(prep$groupModels,do.call,what=hessian_gaussian_kappa_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    h_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * h_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",h_per_group)
}
