# Fast duplication matrix:
duplicationMatrix <- function(n, diag = TRUE){
  mat <- matrix(0,n,n)
  mat[lower.tri(mat,diag=diag)] <- seq_len(sum(lower.tri(mat,diag=diag)))
  mat[upper.tri(mat,diag=diag)] <- t(mat)[upper.tri(mat,diag=diag)]
  inds <- c(mat)
  Matrix::sparseMatrix(i=seq_len(n^2)[inds!=0],j=inds[inds!=0],dims = c(n^2,sum(lower.tri(mat,diag=diag))))
}

# Fast elimination matrix:
eliminationMatrix <- function(n, diag = TRUE){
  mat <- matrix(0,n,n)
  mat[lower.tri(mat,diag=diag)] <- seq_len(sum(lower.tri(mat,diag=diag)))
  inds <- c(mat)
  Matrix::sparseMatrix(j=seq_len(n^2)[inds!=0],i=inds[inds!=0],dims = c(sum(lower.tri(mat,diag=diag)),n^2))
}

# Diagonalization matrix:
E <- function(i,n){
  E <- matrix(0,n,n)
  E[i,i] <- 1
  E
}

diagonalizationMatrix <- function(n){
    as(Reduce(cbind,lapply(1:n,function(x)c(E(x,n)))),"sparseMatrix")
}