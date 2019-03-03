# Fit function per group:
logLikelihood_gaussian_group <- function(S,kappa,means,mu,sigma,...){
  if (any(eigen(kappa)$values < 0)) {
    kappa <- Matrix::nearPD(kappa)$mat
    if (!all(S==0)){
      SK <- Matrix::nearPD(S %*% kappa)$mat  
    } else {
      SK <- S %*% kappa
    }
  } else {
    SK <- S %*% kappa
  }
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