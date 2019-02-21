# # Fit function per group:
# gradient_gaussian_group <- function(S,kappa,means,mu,D,sigma,...){
#   # Mean part:
#   grad_mean <- -2 * kappa %*% (means - mu)
#   
#   # Network part:
#   grad_kappa <- t(D) %*% (Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma))
#   
#   # Combine and return:
#   rbind(grad_mean,grad_kappa)
# }

# Fit function for the precision: -2n* log likelihood
gradient_lnm <- function(x, model){
  
  # Prepare
  prep <- prepare_lnm(x, model)

  # Jacobian:
  Jgauss <- jacobian_gaussian_sigma(prep)
  
  # d_phi_theta:
  d_phi_theta <- d_phi_theta_lnm(prep)
  
  # Model matrix:
  M <- Mmatrix(model@parameters)
  
  # Full Jacobian:
  Jac <- Jgauss %*% d_phi_theta %*% M
  
  # Return:
  return(as.vector(Jac))
}