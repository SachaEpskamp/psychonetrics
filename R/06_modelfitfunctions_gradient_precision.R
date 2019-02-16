# Fit function per group:
gradient_precision_group <- function(S,kappa,means,mu,D,sigma,...){
  # Mean part:
  grad_mean <- -2 * kappa %*% (means - mu)
  
  # Network part:
  grad_kappa <- t(D) %*% (Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma))
  
  # Combine and return:
  rBind(grad_mean,grad_kappa)
}

# Fit function for the precision: -2n* log likelihood
gradient_precision <- function(x, model){
  # Prepare
  prep <- prepare_precision(x, model)

  # Fit function per group:
  gradient_per_group <- lapply(prep$groupModels,do.call,what=gradient_precision_group)
  
  # Bind by row:
  full_gradient <- Reduce("rBind",gradient_per_group)

  grad <- t(prep$groupModels[[1]]$M) %*% full_gradient
  
  # Return:
  return(as.vector(grad))
}