sparseordense <- function(x){

  # Make matrix:
  if (!is.matrix(x)) x <- as.matrix(x)
  
  # Obtain attribute:
  logdet <- attr(x, "logdet")
  
  # Run test for diag - sparse - dense:
  matType <- diag_sparse_dense_cpp(x)
  # 0 = diag, 1 = sparse, 2 = dense
  
  # diagonal:
  if (matType == 0){
    res <- as(x, "dMatrix") 
    
    # Sparse:
  } else if (matType == 1){
    res <- as(x, "dMatrix")
    
    # Else dense:
  } else {
    res <- as.matrix(x)
  }
  # 
  # # First check if diagonal:
  # x2 <- x <- as.matrix(x)
  # diag(x2) <- 0
  # if (all(x2==0)){
  #   res <- as(x, "dMatrix")
  # } else {
  #   # Check if sparse:
  #   if (mean(c(x2[lower.tri(x2,diag=FALSE)],x2[upper.tri(x2,diag=FALSE)]) == 0) > 0.75){
  #     res <- as(x, "dMatrix")
  #   }  else {
  #     res <- as.matrix(x)
  #   }
  # }
  # 
  if (!is.null( logdet)){
    attr(res, "logdet") <- logdet 
  }
  
  
  return(res)
}