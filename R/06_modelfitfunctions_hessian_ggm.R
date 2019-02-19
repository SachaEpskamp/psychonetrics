# # Fit function per group:
# gradient_gaussian_group <- function(S,kappa,means,mu,D,sigma,...){
#   # Mean part:
#   grad_mean <- -2 * kappa %*% (means - mu)
#   
#   # Network part:
#   grad_kappa <- t(D) %*% (Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma))
#   
#   # Combine and return:
#   rBind(grad_mean,grad_kappa)
# }

# Fit function for the precision: -2n* log likelihood
hessian_ggm <- function(x, model){
  
  # Prepare
  prep <- prepare_ggm(x, model)

  # Model matrix:
  M <- Mmatrix(model@parameters)
  
  # Hessian:
  H <- H_ggm(prep)
  
  # Full Hessian:
  Hes <- t(M) %*% H %*% M
  
  # Make symmetric # FIXME: why?
  Hes <- 0.5* (Hes+ t(Hes))
  
  # Return:
  return(as.matrix(Hes))
}