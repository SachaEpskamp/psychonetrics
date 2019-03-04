# Implied model for precision. Requires appropriate model matrices:
implied_gvar_norawts <- function(x){

  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){

    # Implied exogenous covariances:
    exoCov <- x[[g]]$exogenous_sigma

    # Implied stationary distribution (vectorized)
    vecSigma <- x[[g]]$BetaStar %*% x[[g]]$sigmaZetaVec
    Sigma0 <- Matrix(as.vector(vecSigma),nrow = nrow(exoCov), ncol = ncol(exoCov))
    
    # Implied lag-1:
    Sigma1 <- x[[g]]$beta %*% Sigma0
    
    # Full implied sigma:
    sigma <- rbind(
      cbind(exoCov,t(Sigma1)),
      cbind(Sigma1,Sigma0)
    )
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    sigma <- 0.5*(sigma + t(sigma))
    
    # Precision:
    kappa <- solve(sigma)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    kappa <- 0.5*(kappa + t(kappa))
    # Implied variance--covariance:
    # sigma <- as(solve(kappa),"dpoMatrix")
    # sigma <- as(solve(kappa),"dpoMatrix")
    # sigma <- solve(kappa)

    
    return(list(
      kappa = kappa,
      sigma = sigma
      # mu = mu
    )
    )
  })
  
  Result
}
