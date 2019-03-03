# Fit function per group:
maxLikEstimator_Gauss_group <- function(S,kappa,means,mu,sigma,...){
  if (any(eigen(kappa)$values < 0)) {
    kappa <- Matrix::nearPD(kappa)$mat
    # if (!all(S==0)){
    #   SK <- Matrix::nearPD(S %*% kappa)$mat  
    # } else {
    #   SK <- Diagonal(n = nrow(S))
    # }
  } else {
    # if (!all(S==0)){
    #   SK <- S %*% kappa
    # } else {
    #   SK <- Diagonal(n = nrow(S))
    # }
  }

  nvar <- ncol(kappa) 
  res <-  sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu)  - log(det(kappa))
  # res <- sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu) +
    # log(det(kappa)) 
  # res <- sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu)
  as.numeric(res)
}

# Fit function for Gauss ML: -2n* log likelihood
maxLikEstimator_Gauss <- function(x, model){
  # Prepare
  prep <- prepareModel(x, model)

  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=maxLikEstimator_Gauss_group)

  # Sum and return:
  sum(fit_per_group)
}