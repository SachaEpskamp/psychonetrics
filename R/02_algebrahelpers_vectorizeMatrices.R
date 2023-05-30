# Fast duplication matrix:
duplicationMatrix <- function(n, diag = TRUE){
  mat <- matrix(0,n,n)
  mat[lower.tri(mat,diag=diag)] <- seq_len(sum(lower.tri(mat,diag=diag)))
  mat[upper.tri(mat,diag=diag)] <- t(mat)[upper.tri(mat,diag=diag)]
  inds <- c(mat)
  res <- Matrix::sparseMatrix(i=seq_len(n^2)[inds!=0],j=inds[inds!=0],dims = c(n^2,sum(lower.tri(mat,diag=diag))))
  # if (diag){
  #   res <- as(res, "indMatrix")
  # } else {
  #   res <- as(res, "dMatrix")
  # }
  res <- as(res, "dMatrix")
  return(res)
}

# Fast elimination matrix:
eliminationMatrix <- function(n, diag = TRUE){
  mat <- matrix(0,n,n)
  mat[lower.tri(mat,diag=diag)] <- seq_len(sum(lower.tri(mat,diag=diag)))
  inds <- c(mat)
  res <- Matrix::sparseMatrix(j=seq_len(n^2)[inds!=0],i=inds[inds!=0],dims = c(sum(lower.tri(mat,diag=diag)),n^2))
  # res <- as(res, "indMatrix")
  res <- as(res, "dMatrix")
  return(res)
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
  mat <- sparseMatrix(i = seq(1, n^2, length = n),j=seq_len(n),dims = c(n^2,n))
  as(mat, "dMatrix")
}

# Basis matrix:
basisVector <- function(x,n){
  sparseMatrix(i = x,j=1,dims=c(n,1))
}
basisMatrix <- function(n){
  mat <- do.call(rbind,lapply(seq_len(n),function(i,n){
    kronecker(Diagonal(n), basisVector(i,n))
  },n=n))
  as(mat, "dMatrix")
}

# E matrix:
Emat <- function(n,mat = diag(n)){
  mat <- do.call(rbind,lapply(seq_len(n),function(i,n){
    kronecker(Diagonal(n),mat[,i])
  },n=n))
  as(mat, "dMatrix")
}