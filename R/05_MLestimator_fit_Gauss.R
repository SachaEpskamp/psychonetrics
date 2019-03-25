# Fit function per group:
maxLikEstimator_Gauss_group <- function(S,kappa,means,mu,sigma,...){
  
  # nvar <- ncol(kappa) 
  res <-  sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu)  - attr(kappa, "logdet")
  
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