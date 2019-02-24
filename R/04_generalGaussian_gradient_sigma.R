jacobian_gaussian_group_sigmaVersion_meanPart <- function(sigma,mu,means,kappa,...){
  # Mean part:
  # grad_mean <- -2 * t(means - mu) %*% kappa
  grad_mean <- -2 * t(means - mu) %*% kappa
  grad_mean
}

jacobian_gaussian_group_sigmaVersion_sigmaPart <- function(S,means,mu,sigma,D,kappa,...){
  # Network part:
  grad_sigma <- - t(Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma)) %*% (kappa %(x)% kappa) %*% D
  grad_sigma
}


# jacobian function per group:
jacobian_gaussian_group_sigma <- function(...){
  # Mean part:
  grad_mean <- jacobian_gaussian_group_sigmaVersion_meanPart(...)
  
  # Network part:
  grad_sigma <- jacobian_gaussian_group_sigmaVersion_sigmaPart(...)

  # Combine and return:
  cbind(grad_mean,grad_sigma)
}

# Now for all groups:
jacobian_gaussian_sigma <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_gaussian_group_sigma)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
