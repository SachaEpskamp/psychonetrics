# Implied model for precision. Requires appropriate model matrices:
implied_var1 <- function(model,all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }

  if (model@cpp){
    x <- impliedcovstructures_cpp(x,"zeta",type = model@types$zeta, all = all)
  } else {
    x <- impliedcovstructures(x,"zeta",type = model@types$zeta, all = all)
  }
  
  
  
  # For each group:
  nGroup <- length(x)
  for (g in seq_len(nGroup)){
  

    # Some stuff needed now:
    BetaStar <- as.matrix(solve(Diagonal(nrow(x[[g]]$beta)^2) - (x[[g]]$beta %x% x[[g]]$beta)))
    sigmaZetaVec <- Vec(x[[g]]$sigma_zeta)
    
    # Implied exogenous covariances:
    exoCov <- x[[g]]$exo_cholesky %*% t( x[[g]]$exo_cholesky)
    
    # Implied stationary distribution (vectorized)
    vecSigma <- BetaStar %*% sigmaZetaVec
    Sigma0 <- matrix(as.vector(vecSigma),nrow = nrow(exoCov), ncol = ncol(exoCov))
    
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
    if (!sympd_cpp(x[[g]]$sigma)){
      stop("Non positive-definite implied variance-covariance matrix found. Aborting.")
    }
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
    
    # # Kappa, sigma and mu never sparse:
    # x[[g]]$mu <- as.matrix(x[[g]]$mu)
    # x[[g]]$kappa <- as.matrix(x[[g]]$kappa)
    # x[[g]]$sigma <- as.matrix(x[[g]]$sigma)
    browser() 
  }

  
  x
}
