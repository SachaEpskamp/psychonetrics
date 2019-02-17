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

  # Jacobian:
  Hgauss <- hessian_gaussian_kappa(prep)
  
  # d_phi_theta:
  d_phi_theta <- d_phi_theta_ggm(prep)
  
  # Model matrix:
  M <- Mmatrix(model@parameters)
  
  # Full Hessian: # FIXME: This is not true :(
  browser()
  Hes <- t(M) %*% t(d_phi_theta) %*% Hgauss %*% d_phi_theta %*% M
  
  # Return:
  return(as.matrix(Hes))
}