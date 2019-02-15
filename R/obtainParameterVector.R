
parVector <- function(x){
  # Obtain the starting values:
  start <- numeric(max(x@parameters$par))
  for (i in seq_along(start)){
    start[i] <- x@parameters$est[which(x@parameters$par==i)[1]]
  }
  start
}