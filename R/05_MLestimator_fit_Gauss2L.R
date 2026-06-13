# Two-level (random intercept) sufficient-statistics ML estimator for
# multivariate normal data with cluster structure (ml_lvm, estimator = "ML").
#
# Implemented clean-room from the published two-level likelihood
# decomposition (McDonald & Goldstein, 1989; Muthen, 1990):
#
#   -2 l* = (N - J) [ ln|Sigma_W| + tr(Sigma_W^-1 S_PW) ] +
#           sum_s m_s [ ln|Sigma_s| + n_s tr(Sigma_s^-1 A_s) ],
#
# with Sigma_s = Sigma_W + n_s Sigma_B, A_s = cov.d_s +
# (mean.d_s - mu)(mean.d_s - mu)', S_PW the pooled within-cluster covariance
# (ML denominator N - J), and per distinct cluster size n_s: m_s clusters,
# mean.d_s the average of their cluster means and cov.d_s the ML covariance of
# those cluster means. The 2*pi constant is omitted here so that this fit
# function is numerically identical to the existing wide-format FIML objective
# (fimlEstimator_Gauss) at identical parameter values; see
# logLikelihood_gaussian2L_group for the full log-likelihood.
#
# Phase 4 adds a separate MISSING-DATA branch (minustwo_logl_Gauss2L_missing),
# selected when the two-level statistics for a group carry $missing == TRUE.
# It follows the DESIGN (not a code port) of lavaan's
# lav_mvnorm_cluster_missing.R for the special case in which every observed
# variable is both within- and between-level (the ml_lvm case, so lavaan's
# "between-only" block is empty): per missingness pattern P, the within
# precision Sigma_W^(P)^{-1} and log|Sigma_W^(P)| are obtained by a symmetric
# inverse update of Sigma_W^{-1}; per cluster j these are accumulated into
# A_j = sum_{i in j} Sigma_W^(P_i)^{-1} and p_j = sum_i Sigma_W^(P_i)^{-1}
# (y_i - mu)_observed (zero-padded to length p), giving
#
#   -2 l = sum_P f_P ln|Sigma_W^(P)| + sum_j ln|I + Sigma_B A_j|
#          + sum_i (y_i - mu)' Sigma_W^(P_i)^{-1} (y_i - mu)
#          - sum_j p_j' (I + Sigma_B A_j)^{-1} Sigma_B p_j.

# -2 l* for one COMPLETE-data group, WITHOUT the N*p*log(2*pi) constant:
minustwo_logl_Gauss2L_complete <- function(mu, sigma_within, sigma_between, twolevel){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)

  N <- twolevel$N
  J <- twolevel$J

  # Non positive-definite within covariance: return a large penalty, mirroring
  # maxLikEstimator_Gauss_group:
  if (!sympd_cpp(SW)){
    return(1e20)
  }

  # Pooled-within part (zero weight when all clusters have size 1):
  res <- 0
  if (N > J){
    iSW <- solve_symmetric(SW)
    res <- (N - J) * (as.numeric(determinant(SW, logarithm = TRUE)$modulus) + sum(iSW * twolevel$S_PW))
  }

  # Between part, per distinct cluster size:
  sizes <- twolevel$sizes
  for (s in seq_len(nrow(sizes))){
    nj <- sizes$nj[s]
    m  <- sizes$m[s]
    Sj <- SW + nj * SB
    if (!sympd_cpp(Sj)){
      return(1e20)
    }
    iSj <- solve_symmetric(Sj)
    yc <- twolevel$mean_d[[s]] - mu
    A  <- twolevel$cov_d[[s]] + outer(yc, yc)
    res <- res + m * (as.numeric(determinant(Sj, logarithm = TRUE)$modulus) + nj * sum(iSj * A))
  }

  as.numeric(res)
}

# Shared per-cluster quantities for the MISSING-data two-level likelihood.
# Returns the building blocks used by both the fit function and the analytic
# gradient (so the gradient does not recompute them). On a non-positive-definite
# Sigma_W it returns NULL (the caller then applies a large penalty).
twolevel_missing_core <- function(mu, SW, SB, twolevel){
  p <- twolevel$p
  J <- twolevel$J
  mu <- as.vector(mu)

  if (!sympd_cpp(SW)){
    return(NULL)
  }

  SW.inv <- solve_symmetric(SW)
  SW.logdet <- as.numeric(determinant(SW, logarithm = TRUE)$modulus)

  W.logdet <- 0
  PIJ <- matrix(0, twolevel$N, p)           # per-unit Sigma_W^(P)^{-1} (y_i - mu)
  Yc <- sweep(twolevel$Y, 2, mu)            # raw residuals (NA at missing)
  Alist <- rep(list(matrix(0, p, p)), J)     # A_j = sum_i Sigma_W^(P_i)^{-1}
  Wfull_pat <- vector("list", length(twolevel$patterns)) # full p x p W per pattern

  for (ip in seq_along(twolevel$patterns)){
    pp <- twolevel$patterns[[ip]]
    na <- pp$na.idx
    obs <- which(pp$pat)
    # Fully-missing rows carry no information (empty observed set); their full
    # within precision is the zero matrix, contributing nothing:
    if (length(obs) == 0){
      Wfull_pat[[ip]] <- matrix(0, p, p)
      next
    }
    if (length(na) > 0){
      upd <- twolevel_symmetric_inverse_update(SW.inv, na, SW.logdet)
      Wp <- upd$inv
      W.logdet <- W.logdet + upd$logdet * pp$freq
      di <- Yc[pp$rows, obs, drop = FALSE]
      PIJ[pp$rows, obs] <- di %*% Wp
      Wfull <- matrix(0, p, p)
      Wfull[obs, obs] <- Wp
    } else {
      Wp <- SW.inv
      W.logdet <- W.logdet + SW.logdet * pp$freq
      PIJ[pp$rows, ] <- Yc[pp$rows, , drop = FALSE] %*% Wp
      Wfull <- SW.inv
    }
    Wfull_pat[[ip]] <- Wfull
    # Accumulate A_j for every cluster touched by this pattern:
    cl_pat <- pp$clusters
    for (j in unique(cl_pat)){
      Alist[[j]] <- Alist[[j]] + Wfull * sum(cl_pat == j)
    }
  }

  # p_j = sum_{i in j} Sigma_W^(P_i)^{-1}(y_i - mu)  (rows ordered 1..J):
  PJ <- rowsum(PIJ, twolevel$cluster, reorder = TRUE)

  q.yy.a <- sum(PIJ * Yc, na.rm = TRUE)
  q.yy.b <- numeric(J)
  IBZA.logdet <- numeric(J)
  Minv_list <- vector("list", J)
  MinvBp_list <- vector("list", J)
  for (j in seq_len(J)){
    Aj <- Alist[[j]]
    pj <- PJ[j, ]
    M <- SB %*% Aj
    diag(M) <- diag(M) + 1                  # I + Sigma_B A_j
    tmp <- determinant(M, logarithm = TRUE)
    IBZA.logdet[j] <- as.numeric(tmp$modulus) * as.numeric(tmp$sign)
    Minv <- solve(M)
    Minv_list[[j]] <- Minv
    MinvBp_list[[j]] <- drop(Minv %*% (SB %*% pj))   # (I + Sigma_B A_j)^{-1} Sigma_B p_j
    q.yy.b[j] <- sum(pj * MinvBp_list[[j]])
  }

  m2ll <- (q.yy.a - sum(q.yy.b)) + (W.logdet + sum(IBZA.logdet))

  list(
    m2ll = m2ll,
    SW.inv = SW.inv,
    Wfull_pat = Wfull_pat,
    Alist = Alist,
    PJ = PJ,
    Yc = Yc,
    Minv_list = Minv_list,
    MinvBp_list = MinvBp_list
  )
}

# -2 l for one MISSING-data group, WITHOUT the (#observed)*log(2*pi) constant:
minustwo_logl_Gauss2L_missing <- function(mu, sigma_within, sigma_between, twolevel){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  core <- twolevel_missing_core(mu, SW, SB, twolevel)
  if (is.null(core)){
    return(1e20)
  }
  as.numeric(core$m2ll)
}

# Dispatch on complete vs missing data:
minustwo_logl_Gauss2L_noconstant <- function(mu, sigma_within, sigma_between, twolevel){
  if (isTRUE(twolevel$missing)){
    minustwo_logl_Gauss2L_missing(mu, sigma_within, sigma_between, twolevel)
  } else {
    minustwo_logl_Gauss2L_complete(mu, sigma_within, sigma_between, twolevel)
  }
}

# Fit function per group: (-2 l*) / J (no 2*pi constant):
maxLikEstimator_Gauss2L_group <- function(mu, sigma_within, sigma_between, twolevel, ...){
  minustwo_logl_Gauss2L_noconstant(mu = mu, sigma_within = sigma_within,
                                   sigma_between = sigma_between, twolevel = twolevel) / twolevel$J
}

# Fit function for the two-level Gaussian ML estimator. Total fit is the
# nobs-weighted (= clusters-weighted) average of the per-group fits, exactly
# as for the other estimators:
maxLikEstimator_Gauss2L <- function(prep){
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels, do.call, what = maxLikEstimator_Gauss2L_group)
  sum(fit_per_group)
}
