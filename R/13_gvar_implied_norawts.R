# Implied model for precision. Requires appropriate model matrices:
implied_gvar_norawts <- function(model,...){
  x <- formModelMatrices(model)
  
  # For each group:
  nGroup <- length(x)
  for (g in seq_len(nGroup)){
    ### FIXME: TEMPORARY SOLUTION, WILL BE REMOVED!
    IminO <- Diagonal(nrow(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta
    
    if (any(eigen(IminO)$values < 0)){
      warning("I - Omega_zeta is not positive definite, gradient may not be correct.")
      x[[g]]$OmegaStar <- corpcor::pseudoinverse(spectralshift(IminO))
    }  else {
      x[[g]]$OmegaStar <- corpcor::pseudoinverse(IminO)
    }
    x[[g]]$DeltaOmegaStar <- x[[g]]$delta_zeta %*% x[[g]]$OmegaStar 
    x[[g]]$BetaStar <- corpcor::pseudoinverse(Diagonal(nrow(x[[g]]$beta)^2) - (x[[g]]$beta %(x)% x[[g]]$beta))
    x[[g]]$L_betaStar <- model@extramatrices$L %*%  x[[g]]$BetaStar 
    x[[g]]$IkronBeta <- model@extramatrices$In %(x)% x[[g]]$beta
    x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
    x[[g]]$SigmaZeta <- as.matrix(x[[g]]$delta_zeta %*% x[[g]]$OmegaStar %*% x[[g]]$delta_zeta)
    x[[g]]$sigmaZetaVec <- as.vector(x[[g]]$SigmaZeta)
    ###
    
    # Implied exogenous covariances:
    exoCov <- x[[g]]$exo_cholesky %*% t( x[[g]]$exo_cholesky)
    
    # Implied stationary distribution (vectorized)
    vecSigma <- x[[g]]$BetaStar %*% x[[g]]$sigmaZetaVec
    Sigma0 <- Matrix(as.vector(vecSigma),nrow = nrow(exoCov), ncol = ncol(exoCov))
    
    # Implied lag-1:
    Sigma1 <- x[[g]]$beta %*% Sigma0
    
    # Full implied sigma:
    x[[g]]$sigma <- rbind(
      cbind(exoCov,t(Sigma1)),
      cbind(Sigma1,Sigma0)
    )
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
    
    # Precision:
    x[[g]]$kappa <- corpcor::pseudoinverse(x[[g]]$sigma)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
    # Implied variance--covariance:
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- corpcor::pseudoinverse(kappa)
    
  }

  
  x
}
