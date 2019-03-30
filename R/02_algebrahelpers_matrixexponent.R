'%^%' <- function(A,b, tol = 1e-10) {
  if (b==0){
    return(diag(nrow(A)))
  } else {
    res <- Reduce('%*%',rep(list(A),b)) 
    # Tolerance:
    res[abs(res) < tol] <- 0
    as(res, "Matrix")
  }
}
