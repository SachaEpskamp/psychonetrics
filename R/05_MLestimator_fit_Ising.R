# Fit function per group:
maxLikEstimator_Ising_group <- function(omega,tau,beta,delta,squares,means,responses,nobs,Z,...){
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  delta <- as.vector(delta)

  # Compute Z:
  # Z <- computeZ(graph, thresholds, as.numeric(beta), responses, delta = delta)
  #
  # Compute summary statistics:
  # FIXME: Not nice, will make things double
  v1 <- as.vector(means * nobs)
  v2 <- as.matrix(squares)

  # Hamiltonian including the Blume-Capel quadratic term (diag(v2) = sum_obs x_i^2):
  H <-  (
    - sum(thresholds * v1) +
      sum(delta * diag(v2)) -
      sum((graph * v2)[lower.tri(v2,diag=FALSE)])
  )
  
  # Fml
  Fml <- 2 * log(Z) + 2 * beta * H / nobs
  
  
  # Return:
  return(Fml)
}


# Fit function for Ising ML: -2n* log likelihood
maxLikEstimator_Ising <- function(prep){
  
  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=maxLikEstimator_Ising_group)

  # Sum and return:
  sum(fit_per_group)
}