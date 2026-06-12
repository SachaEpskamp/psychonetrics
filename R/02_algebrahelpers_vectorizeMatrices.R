# Session-level cache for the algebra helper matrices (duplication,
# elimination, diagonalization, commutation). These matrices only depend on
# their dimension, are expensive to rebuild for every (sub)model (the
# baseline and saturated submodels rebuild the exact same matrices as the
# target model), and are S4 Matrix objects with copy-on-modify semantics, so
# sharing a single instance across models is safe.
algebraMatrixCacheEnv <- new.env(parent = emptyenv())

cachedAlgebraMatrix <- function(key, build){
  res <- get0(key, envir = algebraMatrixCacheEnv, inherits = FALSE)
  if (is.null(res)){
    res <- build # force the promise
    assign(key, res, envir = algebraMatrixCacheEnv)
  }
  res
}

# Fast duplication matrix:
duplicationMatrix <- function(n, diag = TRUE){
  cachedAlgebraMatrix(paste0("D_", n, "_", diag), {
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
    res
  })
}

# Fast elimination matrix:
eliminationMatrix <- function(n, diag = TRUE){
  cachedAlgebraMatrix(paste0("L_", n, "_", diag), {
    mat <- matrix(0,n,n)
    mat[lower.tri(mat,diag=diag)] <- seq_len(sum(lower.tri(mat,diag=diag)))
    inds <- c(mat)
    res <- Matrix::sparseMatrix(j=seq_len(n^2)[inds!=0],i=inds[inds!=0],dims = c(sum(lower.tri(mat,diag=diag)),n^2))
    # res <- as(res, "indMatrix")
    res <- as(res, "dMatrix")
    res
  })
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
  cachedAlgebraMatrix(paste0("A_", n), {
    mat <- sparseMatrix(i = seq(1, n^2, length = n),j=seq_len(n),dims = c(n^2,n))
    as(mat, "dMatrix")
  })
}

# Fast sparse commutation matrix. Constructs the commutation matrix K(n, m)
# directly in sparse form, equivalent to (but much faster than)
# as(lavaan::lav_matrix_commutation(n, m), "dMatrix"), which builds an
# O((n*m)^2) dense matrix first. The result is coerced to the exact Matrix
# class that the dense coercion used to produce (ddiMatrix when K is the
# identity, dsCMatrix when square/symmetric, dgCMatrix otherwise), so that
# downstream consumers see identical objects:
commutationMatrix <- function(n, m = n){
  cachedAlgebraMatrix(paste0("C_", n, "_", m), {
    nm <- n * m
    if (n == 1L || m == 1L){
      # K(n, 1) and K(1, m) are identity matrices:
      Matrix::Diagonal(nm)
    } else {
      # Row k = i + (j - 1) * n (element (i, j) of the transpose, for an
      # n x m input matrix) picks element j + (i - 1) * m of vec(X):
      res <- Matrix::sparseMatrix(
        i = seq_len(nm),
        j = as.vector(matrix(seq_len(nm), m, n, byrow = TRUE)),
        x = 1,
        dims = c(nm, nm)
      )
      if (n == m){
        # Square commutation matrices are symmetric:
        res <- Matrix::forceSymmetric(res)
      }
      res
    }
  })
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