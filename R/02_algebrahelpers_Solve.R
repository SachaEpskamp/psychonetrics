solve_symmetric <- function(x, logdet = FALSE,approx=TRUE){
  res <- solve_symmetric_cpp(as.matrix(x), logdet, sqrt(.Machine$double.eps),approx=approx)
  inv <- res$inv
  if(logdet){
    attr(inv, "logdet") <- res$logdet
  }
  return(inv)
  # # SYMMETRIC AND TO S4:
  # # x <- as((x + t(x))/2.0, "Matrix")
  # # Force symmetric:
  # x <- as.matrix((x + t(x))/2.0)
  # 
  # # Eigen-values:
  # ev <- eigen(x, symmetric=TRUE, only.values=TRUE)$values
  # 
  # # Check pos-def:
  # if(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0) {
  #   # inv <-  MPinv(x),"Matrix")
  #   inv <-  MPinv(x)
  # 
  #   # Make sparse:
  #   # inv[abs(inv) < sqrt(.Machine$double.eps)] <- 0
  #   # inv <- as(inv, "Matrix")
  #   inv <- inv
  #   
  #   
  #   # Add log determinant:
  #   if (logdet){
  #     attr(inv, "logdet") <- log(.Machine$double.eps)  
  #   }
  #   
  # } else {
  # 
  #       # If pos-def, continue as usual:
  #   # inv <- inv.chol(x, logdet = logdet)
  #   invres <- solve_symmetric_cpp(x, logdet)
  #   inv <- invres$inv
  #   if (logdet){
  #     attr(inv,"logdet") <- invres$logdet
  #   }
  # }
  # 
  # 
  # 
  # # Make Matrix:
  # # inv <- as(inv, "Matrix")
  # # Sparse or dense:
  # 
  # # inv <- sparseordense(inv)
  # 
  # return(inv)
  # # 
  # # # Try spectral shift:
  # # tryres <- try({
  # #   y <- corpcor::pseudoinverse(spectralshift(x))
  # # }, silent = TRUE)
  # # if (!is(tryres, "try-error")){
  # #   return(y)
  # # }
  # # 
  # # stop("Could not compute matrix inverse")
}