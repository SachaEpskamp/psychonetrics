# Implied model for precision. Requires appropriate model matrices:
implied_varcov <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    kappa <- corpcor::pseudoinverse(x[[g]]$sigma)

    # Implied means
    # mu <- x[[g]]$mu
    
    return(list(
      kappa = kappa
      # mu = mu
    )
    )
  })
  
  Result
}