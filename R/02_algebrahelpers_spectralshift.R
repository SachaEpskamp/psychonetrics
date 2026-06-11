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
  
  # Spectral shift and return. Subtract (min eigenvalue - 0.001) from the
  # diagonal only, so the smallest eigenvalue of the result is 0.001 and the
  # off-diagonal elements are untouched. (A previous version had the
  # parentheses misplaced, subtracting min(ev) from the diagonal but adding
  # 0.001 to *every* element, which could leave the matrix exactly singular.)
  x <- x - Diagonal(n=nrow(x)) * (min(Re(eigen(x)$values)) - 0.001)
  return(as.matrix(x))
}