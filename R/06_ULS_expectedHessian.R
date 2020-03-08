# Total:
expected_hessian_ULS_Gaussian <- function(prep){
  
  # Exp hessian per group:
  h_per_group <- lapply(prep$groupModels, do.call, what=ULS_Gauss_exphes_pergroup)
  
  for (i in seq_along(h_per_group)){
    h_per_group[[i]] <- h_per_group[[i]] * (prep$nPerGroup[i]+1)/(prep$nTotal)
  }
  
  do.call(bdiag,h_per_group)
}


# Fit per group:
ULS_Gauss_exphes_pergroup <- function(means,S,mu,sigma,WLS.W,estimator,...){
  if (estimator == "DWLS"){
    WLS.W <- Diagonal(x = diag(WLS.W))
  }
  
  # # model is already prepared!
  # # Dummy for WLS:
  # # observed statistics:
  # obs <-c(as.vector(means),Vech(S))
  # 
  # # implied statistics:
  # imp <-   c(as.vector(mu),Vech(sigma))
  
  # ULS:
  
  # 2 * Diagonal(n = length(imp))
  2 * WLS.W
}