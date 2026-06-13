# Sufficient statistics for the two-level (random intercept) Gaussian ML
# estimator (ml_lvm with estimator = "ML"), computed once from the long-format
# data. Clean-room implementation from the published two-level likelihood
# decomposition (McDonald & Goldstein, 1989; Muthen, 1990).
#
# For COMPLETE data the statistics compress the data into per-distinct-size
# cluster-mean moments (the fast path; see twolevel_sufficient_statistics
# below). When the data contain within-cluster MISSING values they no longer
# compress across clusters: the likelihood is then computed per missing pattern
# and per cluster (Phase 4), so the raw long data plus a per-pattern index are
# stored instead (twolevel_missing_statistics).
#
# The complete-data list ($missing == FALSE) has:
#  - N:      total number of (level-1) units
#  - J:      number of clusters
#  - S_PW:   pooled within-cluster covariance matrix,
#            (1/(N-J)) sum_j sum_i (y_ij - ybar_j)(y_ij - ybar_j)'
#            (zero matrix when all clusters have size 1)
#  - sizes:  data.frame(nj, m): distinct cluster sizes and their counts
#  - mean_d: list (per distinct size s) with the mean of the m_s cluster means
#  - cov_d:  list (per distinct size s) with the ML covariance matrix of the
#            m_s cluster means (divide by m_s; zero matrix when m_s == 1)
#
# The missing-data list ($missing == TRUE) has:
#  - N:      total number of (level-1) units that have at least one observed value
#  - J:      number of clusters
#  - nel:    total number of observed (non-NA) scalar values (for the log(2*pi)
#            constant in the log-likelihood)
#  - Y:      N x p raw data matrix (with NA at the missing entries)
#  - cluster: integer cluster index (1..J) per row of Y
#  - patterns: list, one element per distinct missingness pattern, each with
#       $pat (logical length p, TRUE = observed), $na.idx (missing column
#       indices), $rows (row indices of Y in this pattern), $freq (#rows) and
#       $clusters (cluster index per row of $rows)

# Detect whether two-level data for a group contain within-cluster missingness
# (used to choose the compressed complete-data path or the missing-data path):
twolevel_anyNA <- function(Y){
  anyNA(as.matrix(Y))
}

twolevel_sufficient_statistics <- function(Y, cluster){
  Y <- as.matrix(Y)
  if (anyNA(Y)){
    # Missing data: store the per-pattern / raw structures instead of the
    # compressed cluster-mean moments (Phase 4):
    return(twolevel_missing_statistics(Y, cluster))
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
    missing = FALSE,
    N = N,
    J = J,
    S_PW = S_PW,
    sizes = data.frame(nj = sizes, m = m),
    mean_d = mean_d,
    cov_d = cov_d
  )
}

# Per-pattern / raw structures for the two-level Gaussian ML estimator with
# within-cluster missing data (MCAR/MAR). See the missing-data likelihood in
# 05_MLestimator_fit_Gauss2L.R (minustwo_logl_Gauss2L_missing). Rows that are
# entirely missing carry no information and are dropped.
twolevel_missing_statistics <- function(Y, cluster){
  Y <- as.matrix(Y)
  p <- ncol(Y)
  # Keep the original cluster numbering (1..J) so that the number of clusters
  # equals the number of independent observations (matching the wide-format
  # samplestats). Fully-missing rows are KEPT so the cluster index stays
  # aligned; they carry an empty observed set and contribute nothing to the
  # likelihood (their within precision over the empty set is 0x0):
  cl <- as.integer(factor(cluster))
  J <- max(cl)
  N <- nrow(Y)
  R <- !is.na(Y)

  # Distinct missingness patterns (the empty pattern, if any, is harmless):
  patkey <- apply(R, 1, function(r) paste0(which(r), collapse = ","))
  upat <- unique(patkey)
  patterns <- lapply(upat, function(k){
    rows <- which(patkey == k)
    pat <- R[rows[1], ]
    list(
      pat = pat,
      na.idx = which(!pat),
      rows = rows,
      freq = length(rows),
      clusters = cl[rows]
    )
  })

  list(
    missing = TRUE,
    N = N,
    J = J,
    p = p,
    nel = sum(R),
    Y = Y,
    cluster = cl,
    patterns = patterns
  )
}

# Symmetric inverse update (Schur complement): given S.inv = solve(S) and its
# log-determinant S.logdet, return the inverse and log-determinant of the
# submatrix obtained by removing rows/columns rm.idx from S, WITHOUT inverting
# the submatrix from scratch. Same identity used by lavaan's
# lav_matrix_symmetric_inverse_update; derived from the partitioned-inverse
# (Schur) formula. Returns list(inv, logdet).
twolevel_symmetric_inverse_update <- function(S.inv, rm.idx, S.logdet){
  ndel <- length(rm.idx)
  if (ndel == 0L){
    return(list(inv = S.inv, logdet = S.logdet))
  } else if (ndel == 1L){
    h <- S.inv[rm.idx, rm.idx]
    a <- S.inv[-rm.idx, rm.idx, drop = FALSE] / sqrt(h)
    return(list(inv = S.inv[-rm.idx, -rm.idx, drop = FALSE] - tcrossprod(a),
                logdet = S.logdet + log(h)))
  } else if (ndel < NCOL(S.inv)){
    A <- S.inv[rm.idx, -rm.idx, drop = FALSE]
    H <- S.inv[rm.idx, rm.idx, drop = FALSE]
    return(list(inv = S.inv[-rm.idx, -rm.idx, drop = FALSE] - crossprod(A, solve(H, A)),
                logdet = S.logdet + as.numeric(determinant(H, logarithm = TRUE)$modulus)))
  } else {
    # All rows/cols removed (no observed values) -- should not occur after the
    # empty-row drop above:
    return(list(inv = matrix(0, 0, 0), logdet = 0))
  }
}
