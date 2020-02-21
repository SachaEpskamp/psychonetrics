expected_latent_residual_covs <- function(
  covs,
  lambda
){
  nGroup <- length(covs)
  nLatent <- ncol(lambda)
  nNode <- nrow(lambda)
  
  # Compute the expected latent and residual cov matrices:
  expLatSigma <- lapply(1:nGroup,function(x)matrix(0,nLatent,nLatent))
  expResidSigma <- lapply(1:nGroup,function(x)matrix(0,nNode,nNode))
  
  # For each group:
  for (g in 1:nGroup){
    # Current cov estimate:
    curcov <- as.matrix(covs[[g]])
    
    # Cur loadings:
    curLambda <- lambda[,,g]
    if (!is.matrix(curLambda)){
      curLambda <- as.matrix(curLambda)
    }
    
    # Residual variances, let's start by putting the vars on 1/4 times the observed variances:
    Theta <- diag(diag(curcov)/4)
    
    # Check if this is positive definite:
    ev <- eigen(curcov - Theta)$values
    
    # Shrink until it is positive definite:
    loop <- 0
    repeat{
      ev <- eigen(curcov - Theta)$values
      if (loop == 100){
        # give up...
        
        Theta <- diag(nrow(Theta))
        break
      }
      if (all(ev>0)){
        break
      }
      Theta <- Theta * 0.9
      loop <- loop + 1
    }
    
    # Expected residual sigma:
    expResidSigma[[g]] <- Theta
    
    # This means that the factor-part is expected to be:
    factorPart <- curcov - Theta
    
    # Let's take a pseudoinverse:
    inv <- corpcor::pseudoinverse(kronecker(curLambda,curLambda))
    
    # And obtain psi estimate:
    expLatSigma[[g]] <- matrix(inv %*% as.vector(factorPart),nLatent,nLatent)
  }
  return(list(
    latent = expLatSigma,
    residual = expResidSigma
  ))
}