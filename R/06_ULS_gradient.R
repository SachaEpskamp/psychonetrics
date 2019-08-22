# Fit function for Gauss ML: -2n* log likelihood
ULS_gradient_Gauss <- function(prep){
  # gradient per group:
  g_per_group <- lapply(prep$groupModels, do.call, what=ULS_Gauss_gradient_pergroup)
  
  for (i in seq_along(g_per_group)){
    g_per_group[[i]] <- g_per_group[[i]] * (prep$nPerGroup[i]+1)/(prep$nTotal)
  }

  as.vector(Reduce(cbind,g_per_group))
}


# Fit per group:
ULS_Gauss_gradient_pergroup <- function(means,S,mu,sigma,WLS.V,estimator, meanstructure = TRUE, corinput = FALSE,...){

  
  if (estimator == "DWLS"){
    WLS.V <- Diagonal(x = diag(WLS.V))
  }
  
  # # observed statistics:
  # obs <- c(as.vector(means),Vech(S))
  # 
  # # implied statistics:
  # imp <- c(as.vector(mu),Vech(sigma))
  
  # Include means:
  if (meanstructure){
    obs <- as.vector(means)
    imp <- as.vector(mu)
  } else {
    obs <- numeric(0)
    imp <- numeric(0)
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