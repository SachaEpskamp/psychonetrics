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
# A <- function(i,n){
#   A <- matrix(0,n,n)
#   A[i,i] <- 1
#   A
# }
diagonalizationMatrix <- function(n){
  # browser()
  # seq(1, n^2, length = n)
    # foo <- as(Reduce(cbind,lapply(1:n,function(x)c(A(x,n)))),"sparseMatrix")
  sparseMatrix(i = seq(1, n^2, length = n),j=seq_len(n),dims = c(n^2,n))
}

# Basis matrix:
basisVector <- function(x,n){
  sparseMatrix(i = x,j=1,dims=c(n,1))
}
basisMatrix <- function(n){
  do.call(rbind,lapply(seq_len(n),function(i,n){
    kronecker(Diagonal(n), basisVector(i,n))
  },n=n))
}

# E matrix:
Emat <- function(n,mat = diag(n)){
  do.call(rbind,lapply(seq_len(n),function(i,n){
    kronecker(Diagonal(n),mat[,i])
  },n=n))
}