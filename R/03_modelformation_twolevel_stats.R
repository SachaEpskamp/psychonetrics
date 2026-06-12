# Sufficient statistics for the two-level (random intercept) Gaussian ML
# estimator (ml_lvm with estimator = "ML"), computed once from the long-format
# data. Clean-room implementation from the published two-level likelihood
# decomposition (McDonald & Goldstein, 1989; Muthen, 1990).
#
# Returns a list with:
#  - N:      total number of (level-1) units
#  - J:      number of clusters
#  - S_PW:   pooled within-cluster covariance matrix,
#            (1/(N-J)) sum_j sum_i (y_ij - ybar_j)(y_ij - ybar_j)'
#            (zero matrix when all clusters have size 1)
#  - sizes:  data.frame(nj, m): distinct cluster sizes and their counts
#  - mean_d: list (per distinct size s) with the mean of the m_s cluster means
#  - cov_d:  list (per distinct size s) with the ML covariance matrix of the
#            m_s cluster means (divide by m_s; zero matrix when m_s == 1)
twolevel_sufficient_statistics <- function(Y, cluster){
  Y <- as.matrix(Y)
  if (anyNA(Y)){
    stop("Two-level sufficient statistics require complete data (no NA values).")
  }
  N <- nrow(Y)
  p <- ncol(Y)
  cl <- as.integer(factor(cluster))
  J <- max(cl)

  # Cluster sizes and cluster means:
  ns <- tabulate(cl)
  Ybar <- rowsum(Y, cl) / ns

  # Pooled within-cluster covariance (ML denominator N - J):
  if (N > J){
    S_PW <- crossprod(Y - Ybar[cl,,drop=FALSE]) / (N - J)
  } else {
    S_PW <- matrix(0, p, p)
  }

  # Per distinct cluster size:
  sizes <- sort(unique(ns))
  m <- integer(length(sizes))
  mean_d <- vector("list", length(sizes))
  cov_d <- vector("list", length(sizes))
  for (s in seq_along(sizes)){
    idx <- which(ns == sizes[s])
    m[s] <- length(idx)
    Yb <- Ybar[idx,,drop=FALSE]
    mean_d[[s]] <- colMeans(Yb)
    if (m[s] > 1){
      cov_d[[s]] <- crossprod(sweep(Yb, 2, mean_d[[s]])) / m[s]
    } else {
      cov_d[[s]] <- matrix(0, p, p)
    }
  }

  list(
    N = N,
    J = J,
    S_PW = S_PW,
    sizes = data.frame(nj = sizes, m = m),
    mean_d = mean_d,
    cov_d = cov_d
  )
}
