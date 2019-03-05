# Implied model for precision. Requires appropriate model matrices:
implied_rnm <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    sigma <- x[[g]]$Lambda_BetaStar %*% x[[g]]$sigma_zeta %*% t(x[[g]]$Lambda_BetaStar) + x[[g]]$sigma_epsilon
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    sigma <- 0.5*(sigma + t(sigma))
    
    kappa <- corpcor::pseudoinverse(sigma)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    kappa <- 0.5*(kappa + t(kappa))
    # Implied variance--covariance:
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- corpcor::pseudoinverse(kappa)
    
    # Implied means
    mu <- x[[g]]$tau +  x[[g]]$lambda %*%  x[[g]]$BetaStar  %*% x[[g]]$tau_eta
    
    return(list(
      kappa = kappa,
      sigma = sigma,
      mu = mu
    )
    )
  })
  
  Result
}
