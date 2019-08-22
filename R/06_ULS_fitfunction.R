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
ULS_Gauss_pergroup <- function(means,S,mu,sigma,WLS.V,estimator, meanstructure = TRUE, corinput = FALSE,...){
  

  if (estimator == "DWLS"){
     WLS.V <- Diagonal(x = diag(WLS.V))
  }
  
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
  
  # # observed statistics:
  # obs <- c(as.vector(means),Vech(S))
  # 
  # # implied statistics:
  # imp <- c(as.vector(mu),Vech(sigma))
  
  # ULS:
  as.numeric(t(obs - imp) %*% WLS.V %*% (obs - imp))
}