# Fit function per group:
fit_ggm_group <- function(S,kappa,means,mu,...){
  if (any(eigen(kappa)$values < 0)) return(Inf)
  nvar <- ncol(kappa)
  res <- sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu) - 
    log(det(S %*% kappa)) - nvar
  as.numeric(res)
}

# Fit function for the ggm: -2n* log likelihood
fit_ggm <- function(x, model){
  # Prepare
  prep <- prepare_ggm(x, model)

  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fit_ggm_group)

  # Sum and return:
  sum(fit_per_group)
}