fimlEstimator_Gauss_subgroup <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  # Means:
  mu_p <- mu[obs,]
  dat$n *  maxLikEstimator_Gauss_group(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p)
}

# Fit function per group:
fimlEstimator_Gauss_group <- function(fimldata,fulln,sigma,kappa,mu,...){
  # 
  # browser()
  # Sigma <<- sigma
  # Kappa <<- kappa
  # Mu <<- mu
  # Fimldata <<- fimldata
  # 
  # Subgroup models:
  1/fulln * Reduce("+", lapply(fimldata,fimlEstimator_Gauss_subgroup,sigma=sigma,kappa=kappa,mu=mu))
  
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