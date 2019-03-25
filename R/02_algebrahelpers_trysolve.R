trysolve <- function(x){
warning("'trysolve' should be removed!")
  tryres <- try({
    y <- solve(x)
  }, silent = TRUE)
  if (!is(tryres, "try-error")){
    return(y)
  }

  # Try psuedoinverse:
  tryres <- try({
    y <- corpcor::pseudoinverse(x)
  }, silent = TRUE)
  if (!is(tryres, "try-error")){
    return(y)
  }
  
  # Try spectral shift:
  tryres <- try({
    y <- corpcor::pseudoinverse(spectralshift(x))
  }, silent = TRUE)
  if (!is(tryres, "try-error")){
    return(y)
  }

  stop("Could not compute matrix inverse")
}