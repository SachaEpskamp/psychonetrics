logLikelihood_gaussian_group <- function(data,...){
  if (missing(data)){
    return(logLikelihood_gaussian_group_sumstat(...))
  } else {
    return(logLikelihood_gaussian_group_fiml(data=data,...))
  }
}

# FIML version:
logLikelihood_gaussian_group_fiml <-  function(mu,sigma,data,kappa,...){
  n <- nrow(data)

  # Spectral shift kappa (will do nothing if kappa is pos def):
  kappa <- spectralshift(kappa)
  
  # Current fit:
  curF <- 0

  # usual log determinant:
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
    curF <- curF +  (
      logdetK_p - 
        log((2*pi)^(nrow(kappa_p))) - 
        t(y - mu_p) %*% kappa_p %*% (y - mu_p)
    )
  }
  
  as.numeric((1/n) * curF)
}

# Fit function per group:
logLikelihood_gaussian_group_sumstat <- function(S,kappa,means,mu,sigma,...){
  # if (any(eigen(kappa)$values < 0)) {
  #   kappa <- Matrix::nearPD(kappa)$mat
  #   if (!all(S==0)){
  #     SK <- Matrix::nearPD(S %*% kappa)$mat  
  #   } else {
  #     SK <- S %*% kappa
  #   }
  # } else {
  #   SK <- S %*% kappa
  # }
  kappa <- spectralshift(kappa)
  SK <- S %*% kappa
  nvar <- ncol(kappa)
  res <-  log(det(kappa)) - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  as.numeric(res)
}

# Fit function for Gauss ML: -2n* log likelihood
logLikelihood_gaussian <- function(model){
  # Prepare
  prep <- prepareModel(parVector(model), model)

  # Fit function per group:
  ll_per_Group <- prep$nPerGroup/2 * sapply(prep$groupModels,do.call,what=logLikelihood_gaussian_group)

  # Sum and return:
  sum(ll_per_Group)
}