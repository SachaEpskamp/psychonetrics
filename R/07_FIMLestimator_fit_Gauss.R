fimlEstimator_Gauss_subgroup <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  # Means:
  mu_p <- mu[obs]
  dat$n *  maxLikEstimator_Gauss_group(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p)
}

# Fit function per group:
fimlEstimator_Gauss_group <- function(fimldata,fulln,sigma,kappa,mu,...){

  # # 
  # Sigma <<- sigma
  # Kappa <<- kappa
  # Mu <<- mu
  # Fimldata <<- fimldata
  # # 
  # Subgroup models:
  1/fulln * Reduce("+", lapply(fimldata,fimlEstimator_Gauss_subgroup,sigma=sigma,kappa=kappa,mu=mu))
  
}

# C++ version:
fimlEstimator_Gauss_group_cpp_outer<- function(fimldata,fulln,sigma,kappa,mu,...){
  
  fimlEstimator_Gauss_group_cpp(sigma=as.matrix(sigma),mu=as.matrix(mu),kappa = as.matrix(kappa),fimldata = fimldata,epsilon = .Machine$double.eps, n = fulln)   
  
}

# Full version:
fimlEstimator_Gauss_group_cpp_outer_fullFIML <- function(fimldata,fulln,sigma,kappa,mu,...){
  
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
  
  fimlEstimator_Gauss_group_cpp_fullFIML(sigma=sigma,mu=mu,kappa = kappa,fimldata = fimldata,epsilon = .Machine$double.eps, n = fulln) 
  
}



# Fit function for Gauss ML: -2n* log likelihood
fimlEstimator_Gauss <- function(prep){
  
  # Use C++?
  if (prep$cpp){
    if (prep$fullFIML){
      # Fit function per group:
      fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fimlEstimator_Gauss_group_cpp_outer_fullFIML)      
      
    } else {
      # Fit function per group:
      fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fimlEstimator_Gauss_group_cpp_outer)      
    }
    
  } else {
    
    if (prep$fullFIML){
      stop("Full (rowwise) FIML only supported through C++")
    } else {
      # Fit function per group:
      fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fimlEstimator_Gauss_group)    
    }
    
  }
  
  
  # Sum and return:
  sum(fit_per_group)
}