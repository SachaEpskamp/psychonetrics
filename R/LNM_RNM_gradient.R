# Function to obtain some extra matrices for LNM_RNM (only needed once):
extraMatrices <- function(parList){
  # Form the model:
  
  # Number of variables:
  nVar <- nrow(parList$lambda$cur)
  nLat <- ncol(parList$sigma_eta$cur)
  
  D_sigma <- matrixcalc::duplication.matrix(nVar)
  L_sigma <- matrixcalc::elimination.matrix(nVar)    
  
  if (ncol(parList$lambda$cur) > 1){
    C_lambda <- matrixcalc::commutation.matrix(nrow(parList$lambda$cur),ncol(parList$lambda$cur))    
  } else {
    C_lambda <- diag(nrow(parList$lambda$cur))
  }
  
  if (ncol(parList$kappa_eta$cur) > 1){
    D_kappa_eta <- matrixcalc::duplication.matrix(ncol(parList$kappa_eta$cur))    
  } else {
    D_kappa_eta <- matrix(1,1,1) 
  }
  
  
  return(list(
    D_sigma=D_sigma,
    L_sigma=L_sigma,
    C_lambda=C_lambda,
    D_kappa_eta=D_kappa_eta
  ))
}

# Gradient:
Gradient <- function(x,  parList, S, means,benchmark = FALSE,extraMatrices){
  if (missing(extraMatrices)){
    stop("'extraMatrices' may not be missing (run extraMatrices(...))")
  }
  if(benchmark) t0 <- Sys.time()
  # Form the model:
  mod <- formModel(x, parList)
  if(benchmark) {
    t1 <- Sys.time()
    message("Forming model took: ",round(t1-t0,6))
  }
  
  # Number of variables:
  nVar <- ncol(mod$sigma)
  nLat <- ncol(mod$sigma_eta)
  
  # Some duplication and elimination matrices:
  
  ## Part 1: General
  if(benchmark) t0 <- Sys.time()
  # Fml with respect to means:
  d_Fml_mu <- 2 * t(mod$mu-means) %*% mod$kappa
  
  # Fml with respect to varcovs:
  d_Fml_sigma <- -c(S - (means-mod$mu)%*%t(means-mod$mu)) %*% kronecker(mod$kappa, mod$kappa) + t(c(mod$kappa))
  d_Fml_sigma <- d_Fml_sigma %*% extraMatrices$D_sigma
  
  # Fml with respect to phi:
  d1 <- cbind(d_Fml_mu, d_Fml_sigma)
  if(benchmark) {
    t1 <- Sys.time()
    message("Computing first derivative took: ",round(t1-t0,6))
  }
  ## Part 2: Model specific
  if(benchmark) t0 <- Sys.time()
  # d_mu_mu = identitiy matrix
  d_mu_mu <- diag(nVar)
  
  # to save forming an unnesisary large matrix, I combine part 2 and 3 and only retain the columns of free parameters:
  d_mu_mu <- deleteCols(d_mu_mu,"mu",parList)
  if(benchmark) {
    t1 <- Sys.time()
    message("d_mu_mu took: ",round(t1-t0,6))
  }
  
  # derivative of sigma with respect to lambda:
  if(benchmark) t0 <- Sys.time()
  d_sigma_lambda <- extraMatrices$L_sigma %*% (kronecker(mod$lambda %*% mod$sigma_eta, diag(nVar)) + kronecker(diag(nVar), mod$lambda %*% mod$sigma_eta) %*% extraMatrices$C_lambda)
  
  # to save forming an unnesisary large matrix, I combine part 2 and 3 and only retain the columns of free parameters:
  d_sigma_lambda <- deleteCols(d_sigma_lambda,"lambda",parList)
  if(benchmark) {
    t1 <- Sys.time()
    message("d_sigma_lambda took: ",round(t1-t0,6))
  }
  
  # Derivative of sigma with respect to latent network:
  if(benchmark) t0 <- Sys.time()
  d_sigma_kappa_eta <- - extraMatrices$L_sigma %*% kronecker(mod$lambda, mod$lambda) %*% kronecker(mod$sigma_eta, mod$sigma_eta) %*% extraMatrices$D_kappa_eta
  
  # to save forming an unnesisary large matrix, I combine part 2 and 3 and only retain the columns of free parameters:
  d_sigma_kappa_eta <- deleteCols(d_sigma_kappa_eta,"kappa_eta",parList, symmetrical = TRUE)
  
  if(benchmark) {
    t1 <- Sys.time()
    message("d_sigma_kappa_eta took: ",round(t1-t0,6))
  }
  if(benchmark) t0 <- Sys.time()
  # Derivative of sigma with respect to residual network:
  d_sigma_kappa_epsilon <- - extraMatrices$L_sigma %*% kronecker(mod$sigma_epsilon, mod$sigma_epsilon) %*% extraMatrices$D_sigma
  
  # to save forming an unnesisary large matrix, I combine part 2 and 3 and only retain the columns of free parameters:
  d_sigma_kappa_epsilon <- deleteCols(d_sigma_kappa_epsilon,"kappa_epsilon",parList, symmetrical = TRUE)
  if(benchmark) {
    t1 <- Sys.time()
    message("d_sigma_kappa_epsilon took: ",round(t1-t0,6))
  }
  
  # Total number of parameters:
  nTotalPars <- max(unlist(lapply(parList,"[[","par")))
  
  # Number of means:
  nMeans <- length(mod$mu)
  
  # Number of variances and covariances:
  nVarCovs <- sum(lower.tri(mod$sigma, diag = TRUE))
  
  # Empty full step 2 gradient:
  if(benchmark) t0 <- Sys.time()
  d2 <- matrix(0, nrow = nMeans + nVarCovs, ncol = nTotalPars)
  
  # Fill mu:
  d2[1:nMeans,parList$mu$par[parList$mu$par!=0]] <- d_mu_mu
  
  # Fill lambda:
  d2[nMeans + 1:nVarCovs,parList$lambda$par[parList$lambda$par!=0]] <- d_sigma_lambda
  
  # Fill latent network:
  d2[nMeans + 1:nVarCovs,parList$kappa_eta$par[parList$kappa_eta$par!=0 & lower.tri(parList$kappa_eta$par,diag=TRUE)]] <- d_sigma_kappa_eta
  
  # Fill residual network:
  d2[nMeans + 1:nVarCovs,parList$kappa_epsilon$par[parList$kappa_epsilon$par!=0 & lower.tri(parList$kappa_epsilon$par,diag=TRUE)]] <- d_sigma_kappa_epsilon
  if(benchmark) {
    t1 <- Sys.time()
    message("filling D2 took: ",round(t1-t0,6))
  }
  
  # Final gradient:
  if(benchmark) t0 <- Sys.time()
  G <- d1 %*% d2
  if(benchmark) {
    t1 <- Sys.time()
    message("Final result took: ",round(t1-t0,6))
  }
  
  c(G)
}