solve_symmetric <- function(x, logdet = FALSE){
  # SYMMETRIC AND TO S4:
  x <- as((x + t(x))/2.0, "Matrix")
  
  # Eigen-values:
  ev <- eigen(x, symmetric=TRUE, only.values=TRUE)$values
  
  # Check pos-def:
  if(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0) {
    inv <-  as(MPinv(x),"Matrix")
    
    # Make sparse:
    inv[abs(inv) < sqrt(.Machine$double.eps)] <- 0
    inv <- as(inv, "Matrix")
    
    # Add log determinant:
    if (logdet){
      attr(inv, "logdet") <- log(.Machine$double.eps)  
    }
    
  } else {
    # If pos-def, continue as usual:
    inv <- inv.chol(x, logdet = logdet)
  }
  

  
  # Make Matrix:
  # inv <- as(inv, "Matrix")
  
  return(inv)
  # 
  # # Try spectral shift:
  # tryres <- try({
  #   y <- corpcor::pseudoinverse(spectralshift(x))
  # }, silent = TRUE)
  # if (!is(tryres, "try-error")){
  #   return(y)
  # }
  # 
  # stop("Could not compute matrix inverse")
}