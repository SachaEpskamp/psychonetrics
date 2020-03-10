# From graphicalVAR:
computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  # x <- as.matrix(Matrix::forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  if (ncol(beta) == nrow(beta)+1){
    beta <- beta[,-1,drop=FALSE]
  }
  sigma <- solve_symmetric(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}
