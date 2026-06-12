
# Helper: index of the first parameter-table row for each parameter number
# 1..max(par). match() returns the FIRST match, exactly as the former
# which(par == i)[1] per-parameter loop did, but in a single O(nrow) pass
# instead of O(npar * nrow):
firstParRows <- function(x){
  match(seq_len(max(x@parameters$par)), x@parameters$par)
}

parVector <- function(x){
  # Obtain the starting values:
  x@parameters$est[firstParRows(x)]
}

lowerBound <- function(x){
  x@parameters$minimum[firstParRows(x)]
}

upperBound <- function(x){
  x@parameters$maximum[firstParRows(x)]
}

rowMat <- function(x){
  as.numeric(x@parameters$row[firstParRows(x)])
}
colMat <- function(x){
  as.numeric(x@parameters$col[firstParRows(x)])
}
parMat <- function(x){
  rows <- firstParRows(x)
  if (length(rows) == 0L) return(numeric(0))
  x@parameters$matrix[rows]
}
