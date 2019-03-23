
parVector <- function(x){

  # Obtain the starting values:
  start <- numeric(max(x@parameters$par))
  for (i in seq_along(start)){
    start[i] <- x@parameters$est[which(x@parameters$par==i)[1]]
  }
  start
}

lowerBound <- function(x){
  # Obtain the starting values:
  lower <- numeric(max(x@parameters$par))
  for (i in seq_along(lower)){
    lower[i] <- x@parameters$minimum[which(x@parameters$par==i)[1]]
  }
  lower
}

upperBound <- function(x){
  # Obtain the starting values:
  upper <- numeric(max(x@parameters$par))
  for (i in seq_along(upper)){
    upper[i] <- x@parameters$maximum[which(x@parameters$par==i)[1]]
  }
  upper
}

rowMat <- function(x){
  # Obtain the starting values:
  row <- numeric(max(x@parameters$par))
  for (i in seq_along(row)){
    row[i] <- x@parameters$row[which(x@parameters$par==i)[1]]
  }
  row
}
colMat <- function(x){
  # Obtain the starting values:
  col <- numeric(max(x@parameters$par))
  for (i in seq_along(col)){
    col[i] <- x@parameters$col[which(x@parameters$par==i)[1]]
  }
  col
}
parMat <- function(x){
  # Obtain the starting values:
  mat <- numeric(max(x@parameters$par))
  for (i in seq_along(mat)){
    mat[i] <- x@parameters$matrix[which(x@parameters$par==i)[1]]
  }
  mat
}