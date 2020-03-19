# Fit function per group:
maxLikEstimator_Ising_group <- function(omega,tau,beta,squares,means,responses,nobs,Z,...){
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  
  # Compute Z:
  # Z <- computeZ(graph, thresholds, as.numeric(beta), responses)
  # 
  # Compute summary statistics:
  # FIXME: Not nice, will make things double
  v1 <- as.vector(means * nobs)
  v2 <- as.matrix(squares)

  
  H <-  (
    - sum(thresholds * v1) - 
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