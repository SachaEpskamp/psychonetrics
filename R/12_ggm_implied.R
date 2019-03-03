# Implied model for precision. Requires appropriate model matrices:
implied_ggm <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    sigma <- x[[g]]$delta %*% x[[g]]$IminOinv %*% x[[g]]$delta
    kappa <- solve(sigma)

    # Implied means
    # mu <- x[[g]]$mu
    
    return(list(
      kappa = kappa,
      sigma = sigma
      # mu = mu
    )
    )
  })
  
  Result
}
