# Gradient of the two-level Gaussian ML fit function (see
# 05_MLestimator_fit_Gauss2L.R) with respect to the distribution parameters
# phi = [mu; vech(Sigma_W); vech(Sigma_B)]. Clean-room implementation from the
# same published formulas. With K_s = Sigma_s^-1 - n_s Sigma_s^-1 A_s Sigma_s^-1:
#
#   d(-2 l*)/d mu       = -2 sum_s m_s n_s Sigma_s^-1 (mean.d_s - mu)
#   d(-2 l*)/d Sigma_W  = sum_s m_s K_s + (N - J)(Sigma_W^-1 - Sigma_W^-1 S_PW Sigma_W^-1)
#   d(-2 l*)/d Sigma_B  = sum_s m_s n_s K_s
#
# Matrix derivatives map to vech by doubling the off-diagonal elements and
# keeping the diagonal (standard duplication-matrix convention, matching the
# vech/elimination ordering used elsewhere in psychonetrics).
#
# Phase 4 adds a MISSING-data branch (jacobian_gaussian2L_group_missing),
# derived clean-room from the same scalar -2 l of
# minustwo_logl_Gauss2L_missing by matrix differentials and validated against
# numDeriv. It reuses the per-cluster quantities returned by
# twolevel_missing_core (so the gradient does not recompute the likelihood).

# Helper: matrix derivative -> vech gradient (double off-diagonal entries):
vech_grad_2L <- function(G){
  G2 <- 2 * G
  diag(G2) <- diag(G)
  G2[lower.tri(G2, diag = TRUE)]
}

# Jacobian (1 x npar row vector) per COMPLETE-data group: (1/J) d(-2 l*)/d phi.
jacobian_gaussian2L_group_complete <- function(mu, sigma_within, sigma_between, twolevel, ...){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)
  p <- length(mu)

  N <- twolevel$N
  J <- twolevel$J

  g_mu <- numeric(p)
  G_SW <- matrix(0, p, p)
  G_SB <- matrix(0, p, p)

  sizes <- twolevel$sizes
  for (s in seq_len(nrow(sizes))){
    nj <- sizes$nj[s]
    m  <- sizes$m[s]
    iSj <- solve_symmetric(SW + nj * SB)
    yc <- twolevel$mean_d[[s]] - mu
    A  <- twolevel$cov_d[[s]] + outer(yc, yc)
    K  <- iSj - nj * iSj %*% A %*% iSj
    g_mu <- g_mu - m * 2 * nj * as.numeric(iSj %*% yc)
    G_SW <- G_SW + m * K
    G_SB <- G_SB + m * nj * K
  }

  if (N > J){
    iSW <- solve_symmetric(SW)
    G_SW <- G_SW + (N - J) * (iSW - iSW %*% twolevel$S_PW %*% iSW)
  }

  grad <- c(g_mu, vech_grad_2L(G_SW), vech_grad_2L(G_SB)) / J

  # Return as 1 x npar row matrix (mirrors jacobian_gaussian_group_sigma):
  matrix(grad, nrow = 1)
}

# Jacobian (1 x npar row vector) per MISSING-data group: (1/J) d(-2 l)/d phi.
# Derived from -2 l = q.yy.a - sum_j q.yy.b + sum_P f_P ln|Sigma_W^(P)| +
# sum_j ln|I + Sigma_B A_j|. Using the symmetry of M_j^{-1} Sigma_B (which
# holds because Sigma_B (I + A_j Sigma_B) = (I + Sigma_B A_j) Sigma_B):
#
#   d/d mu     : -2 sum_j (p_j - A_j M_j^{-1} Sigma_B p_j)
#   d/d Sigma_B: sum_j [ Sigma_B^? ... ] (see below; logdet + quadratic terms)
#   d/d Sigma_W: adjoint through the per-pattern within precisions
#                (dSigma_W^(P)^{-1} = -W_P dSigma_W^(P) W_P) accumulating the
#                sensitivities to A_j (via M_j) and to p_j / the per-unit terms.
jacobian_gaussian2L_group_missing <- function(mu, sigma_within, sigma_between, twolevel, ...){
  SW <- as.matrix(sigma_within)
  SB <- as.matrix(sigma_between)
  mu <- as.vector(mu)
  p <- twolevel$p
  J <- twolevel$J

  core <- twolevel_missing_core(mu, SW, SB, twolevel)
  if (is.null(core)){
    # Non-pd within covariance: return a zero gradient (the fit returns the
    # large penalty 1e20; the optimizer steps away from this region):
    return(matrix(0, nrow = 1, ncol = p + p * (p + 1)))
  }

  PJ <- core$PJ
  Alist <- core$Alist
  Minv_list <- core$Minv_list
  MinvBp_list <- core$MinvBp_list

  ## ---- mu gradient ----
  g_mu <- numeric(p)
  for (j in seq_len(J)){
    g_mu <- g_mu - 2 * (PJ[j, ] - drop(Alist[[j]] %*% MinvBp_list[[j]]))
  }

  ## ---- Sigma_B gradient ----
  # Per cluster, -2 l contributes  ln|I + Sigma_B A_j| - p_j' M_j^{-1} Sigma_B p_j.
  #   d ln|I + Sigma_B A_j| / d Sigma_B = (A_j M_j^{-1})' = M_j^{-T} A_j'
  #   d (p_j' M_j^{-1} Sigma_B p_j) / d Sigma_B
  #       = (M_j^{-T} p_j) p_j'  -  (M_j^{-T} p_j) (A_j M_j^{-1} Sigma_B p_j)'
  # and -2 l carries the quadratic with a minus sign.
  G_SB <- matrix(0, p, p)
  for (j in seq_len(J)){
    Minv <- Minv_list[[j]]
    pj <- PJ[j, ]
    MinvBp <- MinvBp_list[[j]]              # M_j^{-1} Sigma_B p_j
    Mtp <- drop(t(Minv) %*% pj)             # M_j^{-T} p_j
    G_logdet <- t(Alist[[j]] %*% Minv)
    G_q <- outer(Mtp, pj) - outer(Mtp, drop(Alist[[j]] %*% MinvBp))
    G_SB <- G_SB + G_logdet - G_q
  }

  ## ---- Sigma_W gradient (adjoint through per-pattern within precisions) ----
  # Sensitivity of -2 l to A_j (holding p_j and the per-unit terms fixed):
  #   GA_j = Sigma_B M_j^{-T} + (Sigma_B M_j^{-T} p_j)(M_j^{-1} Sigma_B p_j)'
  # Sensitivity of -2 l to p_j (holding A_j fixed): Gp_j = -2 M_j^{-1} Sigma_B p_j.
  GA <- vector("list", J)
  Gp <- vector("list", J)
  for (j in seq_len(J)){
    Minv <- Minv_list[[j]]
    pj <- PJ[j, ]
    MinvBp <- MinvBp_list[[j]]
    GA[[j]] <- SB %*% t(Minv) + (SB %*% (t(Minv) %*% pj)) %*% t(MinvBp)
    Gp[[j]] <- -2 * MinvBp
  }

  G_SW <- matrix(0, p, p)
  for (ip in seq_along(twolevel$patterns)){
    pp <- twolevel$patterns[[ip]]
    o <- which(pp$pat)
    # Fully-missing rows contribute nothing to the gradient:
    if (length(o) == 0) next
    Wp <- core$Wfull_pat[[ip]][o, o, drop = FALSE]   # Sigma_W^(P)^{-1}
    di <- core$Yc[pp$rows, o, drop = FALSE]          # observed residuals (freq x |o|)

    # Sensitivity of -2 l to this pattern's within precision Wp (free):
    GWp <- matrix(0, length(o), length(o))
    cl_pat <- pp$clusters
    for (j in unique(cl_pat)){
      GWp <- GWp + sum(cl_pat == j) * GA[[j]][o, o, drop = FALSE]
    }
    # Through PIJ (-> q.yy.a) and p_j (-> q.yy.b): per unit i the residual-side
    # sensitivity is r_i = Yc_i + Gp_{cluster(i)}; PIJ_i = Wp di contributes r_i di':
    Rmat <- di
    for (kk in seq_along(pp$rows)){
      Rmat[kk, ] <- Rmat[kk, ] + Gp[[cl_pat[kk]]][o]
    }
    GWp <- GWp + crossprod(Rmat, di)                 # sum_i r_i di'

    # Map the Wp-sensitivity to Sigma_W^(P): d Wp = -Wp d Sigma_W^(P) Wp, so the
    # Sigma_W^(P)-gradient is -Wp' GWp Wp'; plus the direct per-pattern
    # log-determinant term f_P ln|Sigma_W^(P)| with gradient f_P Sigma_W^(P)^{-T}:
    G_SWoo <- -t(Wp) %*% GWp %*% t(Wp) + pp$freq * t(Wp)
    G_SW[o, o] <- G_SW[o, o] + G_SWoo
  }

  # Symmetrize the matrix derivatives before mapping to vech (the adjoint
  # assembly produces the correct value but in an unsymmetrized form):
  G_SW <- (G_SW + t(G_SW)) / 2
  G_SB <- (G_SB + t(G_SB)) / 2

  grad <- c(g_mu, vech_grad_2L(G_SW), vech_grad_2L(G_SB)) / J
  matrix(grad, nrow = 1)
}

# Dispatch on complete vs missing data:
jacobian_gaussian2L_group_sigma <- function(mu, sigma_within, sigma_between, twolevel, ...){
  if (isTRUE(twolevel$missing)){
    jacobian_gaussian2L_group_missing(mu = mu, sigma_within = sigma_within,
                                      sigma_between = sigma_between, twolevel = twolevel, ...)
  } else {
    jacobian_gaussian2L_group_complete(mu = mu, sigma_within = sigma_within,
                                       sigma_between = sigma_between, twolevel = twolevel, ...)
  }
}

# Now for all groups (mirrors jacobian_gaussian_sigma):
jacobian_gaussian2L_sigma <- function(prep){
  # Jacobian per group:
  g_per_group <- lapply(prep$groupModels, do.call, what = jacobian_gaussian2L_group_sigma)

  # Weight:
  for (i in seq_along(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }

  # Bind by column and return:
  Reduce("cbind", g_per_group)
}
