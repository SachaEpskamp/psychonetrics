simplestructure <- function(x){
  factors <- sort(unique(x))
  nVar <- length(x)
  nFact <- length(factors)
  lambda <- matrix(0,nVar,nFact)

  for (i in seq_along(x)){
    lambda[i,x[i]==factors] <- 1
  }
  colnames(lambda) <- factors
  return(lambda)
}