# Fit function per group:
fit_lnm_group <- function(S,kappa,means,mu,sigma,...){
  if (any(eigen(kappa)$values < 0)) {
    SK <- Matrix::nearPD(S %*% kappa)$mat
    kappa <- Matrix::nearPD(kappa)$mat
  } else {
    SK <- S %*% kappa
  }
  nvar <- ncol(kappa)
  res <- sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu) - 
    log(det(SK)) - nvar
  as.numeric(res)
}

# Fit function for the lnm: -2n* log likelihood
fit_lnm <- function(x, model){
  # Prepare
  prep <- prepare_lnm(x, model)

  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fit_lnm_group)

  # Sum and return:
  sum(fit_per_group)
}