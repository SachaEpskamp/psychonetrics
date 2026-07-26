# Shared helpers for the psychonetrics tinytest suite.
#
# This file is NOT a test file (tinytest only runs files matching "^test");
# it is sourced from the test files themselves.

# Reproducible multivariate normal draws --------------------------------------
#
# MASS::mvrnorm() transforms standard normal draws with the eigendecomposition
# of Sigma. Eigenvectors are only determined up to a sign, and for repeated
# eigenvalues only up to a rotation of the eigenspace, so which decomposition
# LAPACK returns is implementation dependent. Consequently MASS::mvrnorm() does
# NOT produce the same data on different BLAS/LAPACK builds even with a fixed
# seed, which made hard-coded reference values in these tests fail on the CRAN
# ATLAS/MKL/OpenBLAS and macOS-arm64 checks.
#
# The Cholesky factor of a positive definite Sigma is unique, so this generator
# reproduces bit-for-bit (up to floating point rounding) across platforms. It
# consumes the random stream in the same order as MASS::mvrnorm() and mimics its
# drop-to-vector behaviour for n = 1, so it is a drop-in replacement here.
rmvnorm <- function(n, mu, Sigma) {
  p <- length(mu)
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(Sigma)
  X <- X + rep(mu, each = n)
  if (n == 1) drop(X) else X
}
