logLikelihood_gaussian_group <- function(fimldata,...){
  if (missing(fimldata)){
    return(logLikelihood_gaussian_group_sumstat(...))
  } else {
    return(logLikelihood_gaussian_group_fiml(fimldata=fimldata,...))
  }
}

# FIML version:
logLikelihood_gaussian_group_fiml <- function(fimldata,fulln,sigma,kappa,mu,cpp,fullFIML,...){
  
  # print(cpp)
  # Use C++?
  if (cpp){
    if (fullFIML){
      # Fit function per group:
      res <- 1/fulln * logLikelihood_gaussian_subgroup_fiml_cpp_outer_fullFIML(fimldata=fimldata,sigma=sigma,kappa=kappa,mu=mu)
    } else {
      # Fit function per group:
      res <- 1/fulln * logLikelihood_gaussian_subgroup_fiml_cpp_outer(fimldata=fimldata,sigma=sigma,kappa=kappa,mu=mu)
    }
    
  } else {
    if (fullFIML){
      stop("Full (rowwise) FIML only supported through C++")
    } else {
      # Fit function per group:
      res <- 1/fulln * Reduce("+", lapply(fimldata,logLikelihood_gaussian_subgroup_fiml,sigma=sigma,kappa=kappa,mu=mu))    
    }
    
  }
  
  
  
  return(res)

}

# CPP version
logLikelihood_gaussian_subgroup_fiml_cpp_outer <- function(fimldata,sigma,kappa,mu,...){
 
  logLikelihood_gaussian_subgroup_fiml_cpp(sigma=as.matrix(sigma),mu=as.matrix(mu),kappa = as.matrix(kappa),fimldata = fimldata,epsilon = .Machine$double.eps) 
}

# CPP version fullFIML
logLikelihood_gaussian_subgroup_fiml_cpp_outer_fullFIML <- function(fimldata,sigma,kappa,mu,...){
  
  nDat <- length(fimldata)
  
  # Some checks:
  if (!is.list(sigma)){
    sigma <- lapply(seq_len(nDat), function(x) as.matrix(sigma))
  } else {
    if (length(sigma) != nDat){
      stop("Number of 'sigma' matrices must equal number of rows in data.")
    }  
  }
  
  if (!is.list(mu)){
    mu <- lapply(seq_len(nDat), function(x) as.vector(mu))
  } else {
    if (length(mu) != nDat){
      stop("Number of 'mu' vectors must equal number of rows in data.")
    }  
  }
  
  if (!is.list(kappa)){
    kappa <- lapply(seq_len(nDat), function(x) as.matrix(kappa))
  } else {
    if (length(kappa) != nDat){
      stop("Number of 'kappa' matrices must equal number of rows in data.")
    }  
  }
  
  
  logLikelihood_gaussian_subgroup_fiml_cpp_fullFIML(sigma=sigma,mu=mu,kappa = kappa,fimldata = fimldata,epsilon = .Machine$double.eps) 
}


logLikelihood_gaussian_subgroup_fiml <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  

  
  # Means:
  mu_p <- mu[obs]
  dat$n *  logLikelihood_gaussian_group_sumstat(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p)
}


# Fit function per group:
logLikelihood_gaussian_group_sumstat <- function(S,kappa,means,mu,sigma,...){
  SK <- S %*% kappa
  nvar <- ncol(kappa)
  # res <-  log(det(kappa)) - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  res <-  determinant(kappa)$modulus - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  as.numeric(res)
}

# Fit function for Gauss ML: -2n* log likelihood
logLikelihood_gaussian <- function(prep){

  # Add cpp:
  for (i in seq_along(prep$groupModels)){
    prep$groupModels[[i]]$cpp <- prep$cpp
    prep$groupModels[[i]]$fullFIML <- prep$fullFIML
  }

  # Fit function per group:
  ll_per_Group <- prep$nPerGroup/2 * sapply(prep$groupModels,do.call,what=logLikelihood_gaussian_group)

  # Sum and return:
  sum(ll_per_Group)
}