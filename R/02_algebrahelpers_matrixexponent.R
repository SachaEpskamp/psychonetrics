'%^%' <- function(A,b) {
  if (b==0){
    return(diag(nrow(A)))
  } else {
    res <- Reduce('%*%',rep(list(A),b)) 
    # Tolerance:
    res[abs(res) < 1e-3] <- 0
    as(res, "sparseMatrix")
  }
}
