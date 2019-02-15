# Implied model for GGM. Requires appropriate model matrices:
implied_ggm <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    
    # Implied precision:
    kappa <- x[[g]]$kappa
    
    # Implied variance--covariance:
    sigma <- as(solve(kappa),"dpoMatrix")
    # sigma <- as(solve(kappa),"dpoMatrix")
    
    # Implied means
    mu <- x[[g]]$mu
    
    return(list(
      # kappa_z = kappa,
      sigma = sigma
      # mu = mu
    )
    )
  })
  
  Result
}
