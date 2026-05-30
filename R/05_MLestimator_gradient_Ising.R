# jacobian function per group:
jacobian_Ising_group <- function(omega,tau,beta,delta,squares,means,responses,nobs,exp_v1,exp_v2,exp_H,...){
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  delta <- as.vector(delta)
  
  # Compute likelihood table:
  # LikTable <- IsingSampler::IsingLikelihood(graph,thresholds,beta,responses,potential=TRUE)
  # states <- as.matrix(LikTable[,-1])
  # Z <- sum(LikTable$Potential)
  # probabilities <- LikTable$Potential / Z
  # 
  # # Compute expected sum scores:
  # exp_v1 <- colSums(probabilities * states)
  # 
  # # Compute expected products:
  # exp_v2 <- t(probabilities * states) %*% states
  # 
  # 
  # # Compute expected H:
  # expH <- expHcpp(
  #   states,
  #   probabilities,
  #   graph,
  #   thresholds,
  #   nrow(states),
  #   ncol(states)) 
  

  # Compute summary statistics:
  # FIXME: Not nice, will make things double
  v1 <- as.vector(means * nobs)
  v2 <- as.matrix(squares)

  # Hamiltonian (includes the Blume-Capel quadratic term; diag(v2) = sum_obs x_i^2):
  H <-  (
    - sum(thresholds * v1) +
      sum(delta * diag(v2)) -
      sum((graph * v2)[lower.tri(v2,diag=FALSE)])
  )

  # Thresholds gradient:
  threshGrad <- (
    2 * beta * exp_v1 - 2 * beta * v1 / nobs
  )


  # Network: gradient
  graphGrad <- (
    2 * beta * exp_v2 - 2 * beta * v2 / nobs
  )[lower.tri(v2,diag=FALSE)]


  # delta gradient: d Fml / d delta_i = 2 beta ( obs E[x_i^2] - model E[x_i^2] ).
  # delta enters the exponent with a minus sign, hence the opposite sign pattern
  # to the threshold gradient. diag(v2) = sum_obs x_i^2, diag(exp_v2) = E[x_i^2].
  deltaGrad <- (
    2 * beta * diag(v2) / nobs - 2 * beta * diag(exp_v2)
  )

  # beta gradient
  betaGrad <- (
    2*(H/nobs - exp_H)
  )

  # Final gradient (order: tau, omega[lower.tri], delta, beta):
  grad <- Matrix(c(threshGrad, graphGrad, deltaGrad, betaGrad),nrow = 1)
  
  # Return:
  return(grad)

}

# Now for all groups:
jacobian_Ising <- function(prep){
  # model is already prepared!

  # d_phi_theta per group:
  g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_Ising_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
