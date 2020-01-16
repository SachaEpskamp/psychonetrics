# Per group:
expected_hessian_Ising_group <- function(omega,tau,beta,responses,nobs,...){
  
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  
  # Compute likelihood table:
  LikTable <- IsingSampler::IsingLikelihood(graph,thresholds,beta,responses,potential=TRUE)
  states <- as.matrix(LikTable[,-1])
  Z <- sum(LikTable$Potential)
  probabilities <- LikTable$Potential / Z
  
  # Run C++:
  Hes <- expHessianCpp(
    states,
    probabilities,
    graph,
    thresholds,
    beta,
    nrow(states),
    ncol(states))
  
  return(Hes)
  # 
  # # Compute expected sum scores:
  # exp_v1 <- colSums(probabilities * states)
  # 
  # # Compute expected products:
  # exp_v2 <- t(probabilities * states) %*% states
  # 
  # # Compute expected H:
  # expH <- sum(sapply(1:nrow(states),function(i)probabilities[i] * IsingSampler:::H(graph, states[i,], thresholds)))
  # 
  # # Compute block Htt (threholds - thresholds):
  # Htt <- 2 * beta^2 * (exp_v2 - exp_v1 %*% t(exp_v1))
  # 
  # # Compute block Hto (thresholds - omega):
  # Htb <- 
  # 
  # vars <- lapply(1:nrow(states),function(i){
  #   # Probability of this state:
  #   prob <- probabilities[i]
  #   
  #   # Form the vector (v1, v2, H):
  #   vec <- c(states[i,],(states[i,] %*% t(states[i,]))[lower.tri(exp_v2,diag=FALSE)], IsingSampler:::H(graph, states[i,], thresholds)) - 
  #     c(exp_v1, exp_v2[lower.tri(exp_v2,diag=FALSE)], expH)
  #   
  #   # Variance part:
  #   probabilities[i] * (vec %*% t(vec))
  # })
  # 
  # 
  # # Compute the variance:
  # # FIXME: Do this in C++ if you don't want to cry!
  # variance <- as(Reduce("+",vars), "Matrix")
  # 
  # # Return:
  # return(variance)
}

# Total:
expected_hessian_Ising <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_Ising_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",exph_per_group)
}