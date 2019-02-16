# Implied model for precision. Requires appropriate model matrices:
implied_ggm <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    kappa <- x[[g]]$delta %*% (diag(nrow(x[[g]]$omega)) -  x[[g]]$omega) %*% x[[g]]$delta
    
    # Implied variance--covariance:
    # sigma <- as(solve(kappa),"dpoMatrix")
    # sigma <- as(solve(kappa),"dpoMatrix")
    sigma <- solve(kappa)
    
    # Implied means
    mu <- x[[g]]$mu
    
    return(list(
      kappa = kappa,
      sigma = sigma
      # mu = mu
    )
    )
  })
  
  Result
}
