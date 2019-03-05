# Fit function per group:
fimlEstimator_Gauss_group <- function(mu,sigma,data,kappa,...){
  curF <- 0 
  
  n <- nrow(data)

  # Spectral shift kappa (will do nothing if kappa is pos def):
  kappa <- spectralshift(kappa)
  
  # Compute log determinant of kappa:
  logdetK <- log(det(kappa))
  
  # For every subject:
  for (i in seq_len(n)){
    if (!any(is.na(data[i,]))){
      kappa_p <- kappa
      mu_p <- mu
      y <- unlist(data[i,])
      logdetK_p <- logdetK
    } else {
      obs <- unlist(!is.na(data[i,]))
      
      # Skip if nothing to do:
      if (!any(obs)) {
        next
      }
      
      sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
      kappa_p <- corpcor::pseudoinverse(sig_p)

      # Handle possible non positive definiteness:
      kappa_p <- spectralshift(kappa_p)
      
      # Log determinant:
      logdetK_p <- log(det(kappa_p))
      
      # Means:
      mu_p <- mu[obs,]
      
      # raw scores:
      y <- unlist(data[i,obs])
    }
    
    # Add to fit:
    curF <- curF + 1/n * (sum(diag(kappa_p %*% (y - mu_p) %*% t(y - mu_p)))  - logdetK_p)
  }
  
  as.numeric(curF)
}

# Fit function for Gauss ML: -2n* log likelihood
fimlEstimator_Gauss <- function(x, model){
  # Prepare
  prep <- prepareModel(x, model)
  
  
  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fimlEstimator_Gauss_group)
  
  # Sum and return:
  sum(fit_per_group)
}