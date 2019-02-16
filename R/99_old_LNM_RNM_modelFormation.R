formModel <- function(x,parList){
  # Form mu:
  mu <- formMatrix(x,"mu",parList)
  
  # Form lambda:
  lambda <- formMatrix(x,"lambda",parList)
  
  # Form kappa_eta:
  kappa_eta <- formMatrix(x,"kappa_eta",parList)
  
  # Form kappa_epsilon:
  kappa_epsilon <- formMatrix(x,"kappa_epsilon",parList)
  
  # Invert some matrices:
  sigma_eta <- solve(kappa_eta)
  sigma_epsilon <- solve(kappa_epsilon)
  
  # Compute sigma:
  sigma <- lambda %*% sigma_eta %*% t(lambda) + sigma_epsilon
  kappa <- solve(sigma)
  
  # Return everything:
  return(list(
    mu=mu,
    lambda=lambda,
    kappa_eta = kappa_eta,
    sigma_eta = sigma_eta,
    kappa_epsilon = kappa_epsilon,
    sigma_epsilon = sigma_epsilon,
    sigma = sigma,
    kappa = kappa
  ))
}
