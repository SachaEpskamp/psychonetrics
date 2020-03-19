# jacobian function per group:
jacobian_Ising_group <- function(omega,tau,beta,squares,means,responses,nobs,exp_v1,exp_v2,exp_H,...){
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  
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

  # Hamiltionian:  
  H <-  (
    - sum(thresholds * v1) - 
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
  

  # beta gradient
  betaGrad <- (
    2*(H/nobs - exp_H)
  )
  
  # Final gradient:
  grad <- Matrix(c(threshGrad, graphGrad, betaGrad),nrow = 1)
  
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
