# Helpers for the PDC (partial directed correlations) parameterization of
# temporal effects (temporal = "PDC" in var1/panelvar/tsdlvm1;
# temporal_latent / temporal_residual in dlvm1).
#
# The PDC matrix encodes from = row, to = column (the transposed orientation
# relative to beta, matching computePDC() and getmatrix(..., "PDC")):
#
#   PDC[i,j] = beta[j,i] / sqrt(sigma_jj * kappa_ii + beta[j,i]^2),
#
# with sigma the innovation covariance of the temporal process and
# kappa = sigma^{-1}. For fixed innovation structure this map is a smooth
# bijection between beta[j,i] in R and PDC[i,j] in (-1, 1), with explicit
# inverse
#
#   beta[j,i] = PDC[i,j] * sqrt(kappa_ii * sigma_jj) / sqrt(1 - PDC[i,j]^2).
#
# Zeros and signs are preserved, so a (pattern-constrained) PDC model is
# fit-equivalent to the corresponding raw-beta model; equality constraints
# and Wald tests, however, apply directly on the PDC scale.

# Map PDC + innovation covariance -> beta:
PDC_to_beta <- function(PDC, sigma){
  PDC <- as.matrix(PDC)
  sigma <- as.matrix(sigma)
  K <- solve_symmetric(sigma)
  # scale[i,j] = sqrt(kappa_ii * sigma_jj):
  scale <- sqrt(outer(diag(K), diag(sigma)))
  t(PDC * scale / sqrt(pmax(1 - PDC^2, .Machine$double.eps)))
}

# Map beta + innovation covariance -> PDC (same formula as computePDC):
beta_to_PDC <- function(beta, sigma){
  beta <- as.matrix(beta)
  sigma <- as.matrix(sigma)
  K <- solve_symmetric(sigma)
  t(beta / sqrt(outer(diag(sigma), diag(K)) + beta^2))
}

# Reparameterization blocks for the model Jacobian. The Jacobian is
# assembled in raw space (columns for vec(beta) and for the innovation-block
# parameters theta_c); under the PDC parameterization it is post-multiplied
# by the block-triangular map [T X; 0 I]:
#
#   T = d vec(beta) / d vec(PDC)
#     = C %*% diag( sqrt(kappa_ii sigma_jj) * (1 - PDC[i,j]^2)^(-3/2) ),
#   X = d vec(beta) / d theta_c
#     (via d beta[j,i]/d sigma_jj = beta[j,i]/(2 sigma_jj) and
#          d beta[j,i]/d kappa_ii = beta[j,i]/(2 kappa_ii),
#      chained with d vech(sigma)/d theta_c = aug),
#
# so that PDC columns = (beta block) %*% T and the innovation-block columns
# gain (beta block) %*% X. C is the commutation matrix (vec of the transpose)
# and D the duplication matrix of the innovation block:
PDC_reparam <- function(PDC, beta, sigma, aug, D, C){
  PDC <- as.matrix(PDC)
  beta <- as.matrix(beta)
  sigma <- as.matrix(sigma)
  n <- nrow(sigma)
  K <- solve_symmetric(sigma)
  sdiag <- diag(sigma)
  kdiag <- diag(K)

  # T:
  scale <- sqrt(outer(kdiag, sdiag)) # [i,j] = sqrt(kappa_ii sigma_jj)
  Tmat <- C %*% Diagonal(x = as.vector(scale * (1 - PDC^2)^(-3/2)))

  # X:
  Dm <- as.matrix(D)
  # d sigma_jj / d vech(sigma): duplication-matrix rows of the diagonal elements:
  diag_rows <- (seq_len(n) - 1) * n + seq_len(n)
  dsig <- Dm[diag_rows, , drop = FALSE]
  # d kappa_ii / d vech(sigma) = -(K[i,] kron K[i,]) %*% D:
  dkap <- matrix(0, n, ncol(Dm))
  for (i in seq_len(n)){
    dkap[i, ] <- -as.vector(matrix(K[i, ] %x% K[i, ], nrow = 1) %*% Dm)
  }
  X <- matrix(0, n^2, ncol(Dm))
  for (i in seq_len(n)){
    for (j in seq_len(n)){
      b <- beta[j, i]
      if (b != 0){
        # vec index of beta[j,i]:
        X[(i - 1) * n + j, ] <- b / (2 * sdiag[j]) * dsig[j, ] + b / (2 * kdiag[i]) * dkap[i, ]
      }
    }
  }

  list(T = Tmat, X = X %*% as.matrix(aug))
}

# Matrix setup for a PDC parameter matrix. Like matrixsetup_beta, but with
# from = row / to = column orientation ("->"), bounds in (-1, 1), and start
# values obtained by transforming beta start values through beta_to_PDC:
matrixsetup_PDC <- function(
  PDC, # PDC argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "PDC",
  betastart, # optional list of beta start matrices (beta orientation)
  expcov # optional list of expected innovation covariance matrices
){
  # Translate the "diag" shorthand (only autoregressions free):
  if (is.character(PDC) && length(PDC) == 1 && PDC == "diag"){
    PDC <- diag(nNode)
  }

  # Fix PDC:
  PDC <- fixMatrix(PDC, nGroup = nGroup, nrows = nNode, ncols = nNode, equal = equal)

  # For each group, form starting values:
  PDCStart <- PDC
  for (g in seq_len(nGroup)){
    if (!missing(betastart) && !is.null(betastart) && !missing(expcov)){
      sigma <- as.matrix(spectralshift(expcov[[g]]))
      Pst <- beta_to_PDC(as.matrix(betastart[[g]]), sigma)
      # Guard against non-finite transforms of poor starting values:
      Pst[!is.finite(Pst)] <- 0
      # Stay well within the (-1, 1) bounds:
      Pst <- pmax(pmin(Pst, 0.9), -0.9)
      PDCStart[,,g] <- 1*(PDC[,,g]!=0) * Pst
    } else {
      PDCStart[,,g] <- 0.001*(PDC[,,g]!=0)
    }
  }

  # Form the model matrix part:
  list(PDC,
       mat =  name,
       op =  "->",
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       lower = -1,
       upper = 1,
       start = PDCStart,
       sampletable=sampletable
  )
}
