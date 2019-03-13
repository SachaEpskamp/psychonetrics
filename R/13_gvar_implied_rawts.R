# Implied model for precision. Requires appropriate model matrices:
implied_gvar_rawts <- function(x,...){


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
      betak <- x[[g]]$beta%^%l
      as.matrix(betak %*% Sigma0)
    })

    # Now make the Toeplitz block matrix of computational terror:
    sigFull <- as(blockToeplitz(SigmaK),"sparseMatrix")
    
    # And create the massive MU:
    muFull <- Reduce("rbind",lapply(seq_len(nrow(missings)),function(i)x[[g]]$mu))

    # Now cut out all the missing rows:
    # Cut out the rows and cols:
    obsvec <- !as.vector(t(missings))
    muFull <- muFull[obsvec,,drop=FALSE]
    sigFull <- as(sigFull[obsvec,obsvec],"sparseMatrix")
  
    # Precision:


    # This is intractible. let's instead do glasso
    # zeroes <- which(glasso::glasso(as.matrix(sigFull), rho = 0.01)$wi == 0,arr.ind = TRUEwarnings)
    # if (nrow(zeroes) > 0){
    #   kappa <- as(glasso::glasso(as.matrix(sigFull), rho = 0, zero = zeroes)$wi, "sparseMatrix")  
    # } else {
    #   kappa <- as(corpcor::pseudoinverse(sigFull),"sparseMatrix")
    # }
    

    # Let's make kappa artificially sparse:
    kappa <- as(corpcor::pseudoinverse(sigFull),"sparseMatrix")
    kappa[abs(kappa) < 1e-5 & diag(nrow(kappa))!=1] <- 0
    
    
    return(list(
      kappa = as(kappa,"sparseMatrix"),
      sigma = as(sigFull,"sparseMatrix"),
      mu = muFull
    )
    )
  })
  
  Result
}
