# some helper functions:
Vec <- function(x){
  if (class(x)=="matrix"){
    return(matrix(c(x)))
  } else {
    # Assume Matrix package
    return(Matrix(x@x))
  }
}