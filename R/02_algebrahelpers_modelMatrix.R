# Form the model matrix M :
Mmatrix <- function(x){
  rows <- seq_len(nrow(x))[x$par!=0]
  cols <- x$par[x$par!=0]
  dims <- c(nrow(x),max(x$par))
  # Check for diagonal matrix:
  # if (all(rows == cols) && dims[1] == dims[2]){
  #   M <- Diagonal(dims[1]) # Diagonal matrix
  # } else {
    M <- sparseMatrix(i = rows, j = cols, dims = dims)
    M <- as(M, "dMatrix")
  # }
  return(M)
}