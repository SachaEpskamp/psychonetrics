# Implied model for precision. Requires appropriate model matrices:
implied_cholesky <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    # Varcov:
    sigma <- x[[g]]$lowertri %*% t(x[[g]]$lowertri)
    
    # Implied precision:
    kappa <- solve(sigma)

    # Implied means
    # mu <- x[[g]]$mu
    
    return(list(
      sigma = sigma,
      kappa = kappa
      # mu = mu
    )
    )
  })
  
  Result
}
