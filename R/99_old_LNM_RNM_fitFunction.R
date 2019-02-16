# Fit function:
Fml <- function(x, parList, S, means,...){
  # Form the model:
  mod <- formModel(x, parList)
  
  # Number of variables:
  Nvar <- ncol(mod$sigma)
  
  # Fit function:
  c(sum(diag(S %*% mod$kappa)) + 
      t(means-mod$mu) %*% mod$kappa %*% (means-mod$mu) - 
      log(det(S %*% mod$kappa)) - Nvar)
}