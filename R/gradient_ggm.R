# Fit function per group:
gradient_ggm_group <- function(S,kappa,means,mu,D,sigma,...){
  # Mean part:
  grad_mean <- -2 * kappa %*% (means - mu)
  
  # Network part:
  grad_kappa <- t(D) %*% (Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma))
  
  # Combine and return:
  rBind(grad_mean,grad_kappa)
}

# Fit function for the ggm: -2n* log likelihood
gradient_ggm <- function(x, model){
  # Prepare
  prep <- prepare_ggm(x, model)

  # Fit function per group:
  gradient_per_group <- lapply(prep$groupModels,do.call,what=gradient_ggm_group)
  
  # Bind by row:
  full_gradient <- Reduce("rBind",gradient_per_group)
  
  grad <- t(model@fitfunctions$extramatrices$M) %*% full_gradient
  
  # Return:
  return(as.vector(grad))
}