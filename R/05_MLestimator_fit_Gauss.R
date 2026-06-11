# Fit function per group:
maxLikEstimator_Gauss_group <- function(S,kappa,means,mu,sigma,...){
    # nvar <- ncol(kappa)
  # If kappa is not positive definite (e.g. a pseudo-inverse of an indefinite
  # implied covariance matrix), the trace and Mahalanobis terms are unbounded
  # from below and optimizers can diverge to -Inf. Return a large penalty
  # instead, mirroring the C++ twin (maxLikEstimator_Gauss_group_cpp):
  if (!sympd_cpp(as.matrix(kappa))){
    return(1e20)
  }
  # res <-  sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu)  - log(det(kappa))
  res <-  sum(diag(S %*% kappa)) + t(means - mu) %*% kappa %*% (means - mu) - determinant(kappa, TRUE)$modulus

  as.numeric(res)
}

# Fit function for Gauss ML: -2n* log likelihood
maxLikEstimator_Gauss <- function(prep){

    # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=maxLikEstimator_Gauss_group)

  # Sum and return:
  sum(fit_per_group)
}
