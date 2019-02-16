# function to scale from fit function to log likelihood:
loglik_precision <- function(x){
  # Prepare

  n <- sum(x@sample@groups$nobs)
  m <- nrow(x@sample@variables)
  
  # Log det S:
  w <- x@sample@groups$nobs / sum(x@sample@groups$nobs)
  ldS <- sum(w * sapply(x@sample@covs,function(x)log(det(x))))

  -n/2 * (x@objective + m + log((2*pi)^m) + ldS)
}