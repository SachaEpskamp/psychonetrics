# Fit function for Gauss ML: -2n* log likelihood
ULS_gradient_Gauss <- function(prep){
  # gradient per group:
  g_per_group <- lapply(prep$groupModels, do.call, what=ULS_Gauss_gradient_pergroup)
  
  for (i in seq_along(g_per_group)){
    g_per_group[[i]] <- g_per_group[[i]] * (prep$nPerGroup[i]+1)/(prep$nTotal)
  }

  Reduce(cbind,g_per_group)
  # as.vector(Reduce(cbind,g_per_group))
}


# Fit per group:
ULS_Gauss_gradient_pergroup <- function(means,S,mu,tau,thresholds,sigma,WLS.V,estimator, meanstructure = TRUE, corinput = FALSE,...){

  
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

  # ULS:
  -2 * t(obs - imp) %*% WLS.V
}