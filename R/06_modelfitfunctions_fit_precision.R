# Fit function per group:
fit_precision_group <- function(S,kappa,means,mu,...){
  nvar <- ncol(kappa)
  res <- sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu) - 
    log(det(S %*% kappa)) - nvar
  as.numeric(res)
}

# Fit function for the precision: -2n* log likelihood
fit_precision <- function(x, model){
  # Prepare
  prep <- prepare_precision(x, model)

  # Fit function per group:
  fit_per_group <- (prep$nPerGroup / prep$nTotal) * sapply(prep$groupModels,do.call,what=fit_precision_group)

  # Sum and return:
  sum(fit_per_group)
}