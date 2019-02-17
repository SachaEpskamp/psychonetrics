# jacobian function per group:
jacobian_gaussian_group_kappa <- function(S,kappa,means,mu,D,sigma,...){
  # Mean part:
  grad_mean <- -2 * t(means - mu) %*% kappa
  
  # Network part:
  grad_kappa <- t(Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma)) %*% D

  # Combine and return:
  cbind(grad_mean,grad_kappa)
}

# Now for all groups:
jacobian_gaussian_kappa <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_gaussian_group_kappa)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
