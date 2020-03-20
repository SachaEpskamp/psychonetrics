sparseordense <- function(x){
  logdet <- attr(x, "logdet")
  
  # First check if diagonal:
  x2 <- x <- as.matrix(x)
  diag(x2) <- 0
  if (all(x2==0)){
    res <- as(x, "dgCMatrix")
  } else {
    # Check if sparse:
    if (mean(c(x2[lower.tri(x2,diag=FALSE)],x2[upper.tri(x2,diag=FALSE)]) == 0) > 0.75){
      res <- as(x, "dgCMatrix")
    }  else {
      res <- as.matrix(x)
    }
  }
  
  if (!is.null( logdet)){
    attr(res, "logdet") <- logdet 
  }
  
  
  return(res)
}