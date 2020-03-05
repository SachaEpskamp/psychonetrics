spectralshift <- function(x){
# warning("'spectralshift' should be removed!")
  # # Make matrix:
  # if (!is.matrix(x)){
  #   x <- as.matrix(x)
  # }
  # 
  # # Make symmetric:
  # if (!all(x == t(x))){
  #   x <- 0.5 * (x + t(x))
  # }
  # 
  # If anything NA, just return an identity matrix:
  if (any(!is.finite(x))){
    return(Diagonal(n=nrow(x)))
  }
  
  # If all eigenvalues are good, stop:
  if (!any(Re(eigen(x)$values) < 0)){
    return(as.matrix(x))
  }
  
  # Spectral shift and return:
  x <- x - (Diagonal(n=nrow(x)) * (min(Re(eigen(x)$values)))-0.001)
  return(as.matrix(x))
}