# Maximum likelihood covariance estimate:
covML <- function(x,...){
  (nrow(x)-1)/(nrow(x)) * cov(x, ...)
}

# transform between ML and unbiased estimators:
covMLtoUB <- function(x,n){
  n/(n-1) * cov(x, ...)
}

covUBtoML <- function(x,n){
  (n-1)/n * x
}