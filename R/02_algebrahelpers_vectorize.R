# some helper functions:
Vec <- function(x){
  if (class(x)=="matrix"){
    return(matrix(c(x)))
  } else {
    # Assume Matrix package
    if (is(x,"sparseMatrix")){
      return(c(as.matrix(x)))
    } else {
      return(Matrix(x@x))  
    }
  }
}

Vech <- function(x,diag=TRUE){
  if (class(x)!="matrix"){
    x <- as.matrix(x)
  }
  return(x[lower.tri(x,diag=diag)])
}