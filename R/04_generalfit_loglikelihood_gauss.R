logLikelihood_gaussian_group <- function(fimldata,...){
  if (missing(fimldata)){
    return(logLikelihood_gaussian_group_sumstat(...))
  } else {
    return(logLikelihood_gaussian_group_fiml(fimldata=fimldata,...))
  }
}

# FIML version:
logLikelihood_gaussian_group_fiml <- function(fimldata,fulln,sigma,kappa,mu,...){
  # Subgroup models:
  1/fulln * Reduce("+", lapply(fimldata,logLikelihood_gaussian_subgroup_fiml,sigma=sigma,kappa=kappa,mu=mu))
}

logLikelihood_gaussian_subgroup_fiml <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  

  
  # Means:
  mu_p <- mu[obs,]
  dat$n *  logLikelihood_gaussian_group_sumstat(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p)
}


# Fit function per group:
logLikelihood_gaussian_group_sumstat <- function(S,kappa,means,mu,sigma,...){

  SK <- S %*% kappa
  nvar <- ncol(kappa)
  res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
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