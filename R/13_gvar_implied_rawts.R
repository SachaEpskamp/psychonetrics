# Implied model for precision. Requires appropriate model matrices:
implied_gvar_rawts <- function(x){


  # For each group:
  nGroup <- length(x)
  Result <- lapply(seq_len(nGroup), function(g){
    missings <- x[[g]]$mis
    
    # Implied stationary distribution (vectorized)
    vecSigma <- x[[g]]$BetaStar %*% x[[g]]$sigmaZetaVec
    Sigma0 <- Matrix(as.vector(vecSigma),nrow = nrow(x[[g]]$beta), ncol = nrow(x[[g]]$beta))
    
    # How many lag matrices do I need?
    lags <- seq_len(nrow(missings))
    
    # Implied lag-k:
    SigmaK <- lapply(lags-1, function(l){
      as.matrix((x[[g]]$beta%^%l) %*% Sigma0)
    })
    
    # Now make the Toeplitz block matrix of computational terror:
    sigFull <- as(blockToeplitz(SigmaK),"dsyMatrix")
    
    # And create the massive MU:
    muFull <- Reduce("rbind",lapply(seq_len(nrow(missings)),function(i)x[[g]]$mu))

    # Now cut out all the missing rows:
    # Cut out the rows and cols:
    obsvec <- !as.vector(t(missings))
    muFull <- muFull[obsvec,,drop=FALSE]
    sigFull <- sigFull[obsvec,obsvec]
    
    # Precision:

    # This is intractible. let's instead do glasso
    # kappa <- solve(sigFull)
    kappa <- as(glasso::glasso(as.matrix(sigFull), rho = 0.1)$wi, "sparseMatrix")
    
    
    return(list(
      kappa = kappa,
      sigma = as(sigFull,"sparseMatrix"),
      mu = muFull
    )
    )
  })
  
  Result
}
