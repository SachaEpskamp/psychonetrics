# some helper functions:
Vec <- function(x){
  if (is.null(dim(x))){
    return(Matrix(x,1,1))
  }
  
  if (class(x)=="matrix"){
    return(matrix(c(x)))
  } else {
    # Assume Matrix package
    if (is(x,"sparseMatrix")){
      as(c(as.matrix(x)), "Matrix") # FIXME: This should be better...
      # if (length(x@x) == 0){
      #   return(Matrix(0, nrow=nrow(x)*ncol(x),ncol=1))
      # } else {
      #   as(c(as.matrix(x)), "Matrix") # FIXME: This should be better...
      # }
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