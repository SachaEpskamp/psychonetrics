# Implied model for precision. Requires appropriate model matrices:
implied_lnm <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    sigma <- x[[g]]$lambda %*% x[[g]]$sigma_eta %*% t(x[[g]]$lambda ) + x[[g]]$sigma_epsilon
    kappa <- corpcor::pseudoinverse(sigma)
    # Implied variance--covariance:
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- corpcor::pseudoinverse(kappa)
    
    # Implied means
    mu <- x[[g]]$tau +  x[[g]]$lambda %*%  x[[g]]$mu_eta
    
    return(list(
      kappa = kappa,
      sigma = sigma,
      mu = mu
    )
    )
  })
  
  Result
}
