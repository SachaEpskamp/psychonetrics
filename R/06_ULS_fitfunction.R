ULSestimator <- function(x, model){
  # What distribution?
  distribution <- model@distribution
  
  # Function to use:
  distFun <- switch(distribution,
                    "Gaussian" = ULS_Gauss)
  
  # Run and return:
  distFun(x, model)
}


# Fit function for Gauss ML: -2n* log likelihood
ULS_Gauss <- function(x, model){

  # Prepare
  prep <- prepareModel(x, model)

  # Fit per group:
  fit_per_group <- (prep$nPerGroup+1)/(prep$nTotal) * sapply(prep$groupModels, do.call, what=ULS_Gauss_pergroup)

  sum(fit_per_group)
}

# Fit per group:
ULS_Gauss_pergroup <- function(means,S,tau,mu,sigma,WLS.V,estimator,thresholds, meanstructure = TRUE, corinput = FALSE,...){

  if (estimator == "DWLS"){
     WLS.V <- Diagonal(x = diag(WLS.V))
  }
  
  # Include means:
  # FIXME: old code I am not messing with now
  if (missing(tau) || all(is.na(as.matrix(tau)))){
    if (meanstructure){
      obs <- as.vector(means)
      imp <- as.vector(mu)
    } else {
      obs <- numeric(0)
      imp <- numeric(0)
    }
    
  } else {
    obs <- numeric(0)
    imp <- numeric(0)
    
    # Mix means and thresholds:
    for (i in seq_along(thresholds)){
      if (!is.na(means[i])){
        obs <- c(obs, means[i])
        imp <- c(imp, mu[i])
      } else {
        obs <- c(obs, thresholds[[i]])
        imp <- c(imp, tau[seq_len(length(thresholds[[i]])),i])
      }
    }
    
  }

  
  # Add variances (if needed) and covariances:
  if (corinput){
    obs <- c(obs, Vech(S, diag = FALSE))
    imp <- c(imp, Vech(sigma, diag = FALSE))
  } else {
    obs <- c(obs, Vech(S))
    imp <- c(imp, Vech(sigma))
  }
  

  
  # # observed statistics:
  # obs <- c(as.vector(means),Vech(S))
  # 
  # # implied statistics:
  # imp <- c(as.vector(mu),Vech(sigma))
  
  # ULS:
  as.numeric(t(obs - imp) %*% WLS.V %*% (obs - imp))
}