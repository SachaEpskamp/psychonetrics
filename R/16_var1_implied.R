# Implied model for precision. Requires appropriate model matrices:
implied_var1 <- function(model,all = FALSE){
  x <- formModelMatrices(model)

  x <- impliedcovstructures(x,"zeta",type = model@types$zeta, all = all)
  
  # For each group:
  nGroup <- length(x)
  for (g in seq_len(nGroup)){
  

    # Some stuff needed now:
    BetaStar <- as(solve(Diagonal(nrow(x[[g]]$beta)^2) - (x[[g]]$beta %x% x[[g]]$beta)),"Matrix")
    sigmaZetaVec <- Vec(x[[g]]$sigma_zeta)
    
    # Implied exogenous covariances:
    exoCov <- x[[g]]$exo_cholesky %*% t( x[[g]]$exo_cholesky)
    
    # Implied stationary distribution (vectorized)
    vecSigma <- BetaStar %*% sigmaZetaVec
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
    x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    # x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
    # Implied variance--covariance:
    # Extra matrices needed in optimization:
    if (!all){
      x[[g]]$L_betaStar <- model@extramatrices$L %*%  BetaStar 
      x[[g]]$IkronBeta <- model@extramatrices$In %x% x[[g]]$beta
      # x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
      x[[g]]$sigmaZetaVec <- sigmaZetaVec
      x[[g]]$BetaStar <- BetaStar
    } else {
      # Add PDC:
      x[[g]]$PDC <- computePDC(x[[g]]$beta,x[[g]]$kappa_zeta)
    }
    
  }

  
  x
}
