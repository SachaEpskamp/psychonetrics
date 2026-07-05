# Satorra-Bentler-family SCALED chi-square DIFFERENCE tests for nested models.
#
# These power the scaled.test columns of compare() when ALL compared models
# carry a robust scaling factor (estimator = MLM/MLMV/MLMVS/MLR or
# DWLS/WLS/ULS with a stored asymptotic Gamma). For an adjacent nested pair
# (M0 nested in M1, df r0 > r1) with UNscaled statistics T0, T1 and scaling
# factors c0, c1, and m = r0 - r1:
#
#   satorra.bentler.2001 : cd = (r0 c0 - r1 c1)/m ; T = (T0 - T1)/cd
#   satorra.bentler.2010 : as 2001 but with c1 replaced by c10, the scaling
#                          factor of M1's model evaluated at M0's estimates
#                          (strictly positive; Satorra & Bentler 2010)
#   satorra.2000         : exact trace-based test (mean-scaled or scaled-shifted),
#                          required for the scaled-shifted estimators (MLMV/WLSMV)
#
# The math is ported from lavaan's lavTestLRT() and verified against it (see
# extra/fable_audit_scripts/pn_plan_robust/04_diff_tests_verify.R).
# References: Satorra & Bentler (2001); Satorra & Bentler (2010); Satorra (2000);
#   Asparouhov & Muthen (2010).


# Per-parameter-row key used to match parameter-table rows across two models
# (same matrix element in the same group). Matched on matrix/row/col/group_id
# exactly as the task specifies.
scaled_diff_parkey <- function(pt){
  paste(pt$matrix, pt$row, pt$col, pt$group_id, sep = "|")
}

# Basis for the orthogonal complement of the column space of a matrix H
# (mirrors lavaan's lav_matrix_orthogonal_complement): the trailing columns of
# the complete Q factor of a QR decomposition. Returns an (nrow(H) x (nrow-rank))
# matrix whose columns are orthonormal and orthogonal to every column of H.
scaled_diff_orth_complement <- function(H){
  QR <- qr(H)
  Q <- qr.Q(QR, complete = TRUE)
  r <- QR$rank
  if (r >= ncol(Q)) return(matrix(0, nrow(Q), 0))
  Q[, (r + 1):ncol(Q), drop = FALSE]
}

# Constraint Jacobian A (m x npar1) of the restriction that nests M0 in M1,
# expressed in M1's free-parameter space (Satorra 2000 "delta" method, reference
# = H1; see lavaan:::lav_test_diff_A). Built from the two models' full sample-
# statistic Jacobians Delta1, Delta0 (means then vech(Sigma) per group, identical
# ordering): H = ginv(Delta1) Delta0 expresses M0's parameter directions inside
# M1's parameter space, and A spans the orthogonal complement of col(H). Returns
# NULL if the Jacobians cannot be built or have incompatible row counts.
scaled_diff_A <- function(m1, m0){
  Delta1 <- tryCatch(build_Delta_full(m1), error = function(e) NULL)
  Delta0 <- tryCatch(build_Delta_full(m0), error = function(e) NULL)
  if (is.null(Delta1) || is.null(Delta0)) return(NULL)
  if (nrow(Delta1) != nrow(Delta0)) return(NULL)
  H <- tryCatch(MASS::ginv(Delta1) %*% Delta0, error = function(e) NULL)
  if (is.null(H)) return(NULL)
  t(scaled_diff_orth_complement(H))
}

# Estimator-agnostic per-group building blocks for a SINGLE fitted model, in the
# group-weighting convention used by compute_scaled_test_core():
#   Delta : list of group model Jacobians (means then vech(Sigma))
#   Wt    : list of WEIGHT blocks   Wtilde_g = f_g * W_g
#           (W_g = V_g for ML, W_g for WLS, diag(W_g) for DWLS, I for ULS)
#   Gt    : list of GAMMA blocks    Gtilde_g = Gamma_g / f_g
#   bread : the matrix whose inverse projects out the fitted directions
#           (expected Fisher information E for ML; DtWD for least squares)
# Returns NULL when the model is neither robust ML nor a least-squares fit with a
# stored Gamma. Reuses ml_robust_components() / wls_sandwich_components().
scaled_diff_components <- function(x){
  if (x@estimator %in% c("WLS", "DWLS", "ULS")){
    comp <- wls_sandwich_components(x)
    if (is.null(comp)) return(NULL)
    return(list(Delta = comp$Delta, Wt = comp$W, Gt = comp$Gamma,
                bread = comp$DtWD, nGroups = comp$nGroups))
  }
  if (x@estimator == "ML"){
    comp <- ml_robust_components(x)
    if (is.null(comp)) return(NULL)
    fg <- comp$fg; nG <- comp$nGroups
    Wt <- vector("list", nG); Gt <- vector("list", nG)
    for (g in seq_len(nG)){
      Wt[[g]] <- fg[g] * comp$V[[g]]
      Gt[[g]] <- comp$Gamma[[g]] / fg[g]
    }
    return(list(Delta = comp$Delta, Wt = Wt, Gt = Gt,
                bread = comp$E, nGroups = nG))
  }
  # FIML (MLR + missing data) has no fixed-form asymptotic Gamma: the exact
  # Satorra-2000 trace test is not available; the 2001 test (scaling factors
  # only) is still used by the caller.
  NULL
}

# Satorra-Bentler (2010) scaling factor c10: the (mean-adjusted) scaling factor
# of M1's MODEL structure evaluated at M0's estimates, WITHOUT re-optimization
# (Satorra & Bentler 2010; lavaan:::lav_test_diff_m10). M0's estimates are
# injected into M1's parameter table by matching matrix/row/col/group_id; the
# model-implied moments and the expected information are then recomputed at that
# point and the Satorra-Bentler scaling factor is read off. Returns NULL on any
# failure (the caller then falls back to the 2001 test for that pair).
scaled_diff_c10 <- function(m1, m0){
  if (!(m1@estimator %in% c("ML", "WLS", "DWLS", "ULS"))) return(NULL)
  key1 <- scaled_diff_parkey(m1@parameters)
  key0 <- scaled_diff_parkey(m0@parameters)
  match0 <- match(key1, key0)
  ok <- !is.na(match0)
  if (!any(ok)) return(NULL)

  M10 <- m1
  pt <- M10@parameters
  pt$est[ok] <- m0@parameters$est[match0[ok]]
  M10@parameters <- pt
  # Recompute the implied moments at the injected estimates (no optimization):
  M10 <- tryCatch(updateModel(parVector(M10), M10, updateMatrices = TRUE),
                  error = function(e) NULL)
  if (is.null(M10)) return(NULL)
  # Force the expected Fisher information to be recomputed at the injected
  # moments (the stored E belongs to M1's own solution):
  M10@information <- matrix(numeric(0), 0, 0)

  if (M10@estimator %in% c("WLS", "DWLS", "ULS")){
    res <- tryCatch(compute_wlsmv_correction(M10), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    return(res$scaling.factor.sb)
  }
  res <- tryCatch(compute_ml_scaled_test(M10, chisq = m1@fitmeasures$chisq),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  res$scaling.factor.sb
}

# Satorra (2000) exact difference-test traces trace(U Gamma) and trace((U Gamma)^2)
# for the nested pair (M0 in M1), assembled on the GLOBAL block-diagonal stack of
# M1 (same convention as compute_scaled_test_core, so multi-group cross-group
# constraints are handled correctly). With per-group weights already applied
# (Wt_g = f_g W_g, Gt_g = Gamma_g / f_g):
#   P   = bread (M1's expected information / DtWD)
#   Pi  = stacked M1 Jacobian (means then vech(Sigma))
#   A   = constraint Jacobian (scaled_diff_A)
#   paapaap = P^-1 A' (A P^-1 A')^- A P^-1
#   U   = V Pi paapaap Pi' V   (V = block-diagonal Wt)
#   trace(U Gamma), trace((U Gamma)^2) on the block-diagonal G = bdiag(Gt)
# Returns list(trUG, trUG2) or NULL.
scaled_diff_satorra2000_traces <- function(m1, A){
  comp <- scaled_diff_components(m1)
  if (is.null(comp)) return(NULL)
  P <- comp$bread
  P.inv <- tryCatch(solve(P), error = function(e) NULL)
  if (is.null(P.inv) || any(!is.finite(P.inv))) return(NULL)

  bdiag_dense <- function(L){
    n <- sum(vapply(L, nrow, integer(1)))
    M <- matrix(0, n, n); o <- 0L
    for (B in L){ nb <- nrow(B); M[o + seq_len(nb), o + seq_len(nb)] <- B; o <- o + nb }
    M
  }
  V_all <- bdiag_dense(comp$Wt)
  G_all <- bdiag_dense(comp$Gt)
  Pi <- do.call(rbind, comp$Delta)

  APA <- A %*% P.inv %*% t(A)
  APA.inv <- tryCatch(MASS::ginv(APA), error = function(e) NULL)
  if (is.null(APA.inv)) return(NULL)
  paapaap <- P.inv %*% t(A) %*% APA.inv %*% A %*% P.inv

  U <- V_all %*% Pi %*% paapaap %*% t(Pi) %*% V_all
  UG <- U %*% G_all
  trUG <- sum(diag(UG))
  trUG2 <- sum(UG * t(UG))
  if (!is.finite(trUG) || !is.finite(trUG2)) return(NULL)
  list(trUG = trUG, trUG2 = trUG2)
}

# Scaled chi-square DIFFERENCE statistic for ONE adjacent nested pair.
#   m1, m0 : the two fitted psychonetrics models (m0 nested in m1; r0 > r1)
#   method : "satorra.bentler.2001" | "satorra.bentler.2010" | "satorra.2000"
#   scaled.shifted : (satorra.2000 only) scaled-and-shifted vs mean-scaled
# Returns list(stat, df, pvalue) -- stat = NA (with a warning) when the test
# cannot be computed (cd <= 0, m == 0, or a building block is unavailable).
scaled_diff_test_pair <- function(m1, m0, method = "satorra.bentler.2001",
                                  scaled.shifted = TRUE){
  T1 <- m1@fitmeasures$chisq; r1 <- m1@fitmeasures$df
  T0 <- m0@fitmeasures$chisq; r0 <- m0@fitmeasures$df
  m <- r0 - r1

  # The Satorra-Bentler 2001/2010 mean-scaling cd = (r0 c0 - r1 c1)/m requires
  # the MEAN-ADJUSTED scaling factor c = trace(UG)/df. For MLM / MLR / WLSM that
  # equals chisq.scaling.factor, but for the scaled-shifted (MLMV / WLSMV) and
  # mean-and-variance adjusted (MLMVS) tests chisq.scaling.factor is a DIFFERENT
  # quantity (1/a resp. trUG2/trUG). Prefer the separately-stored
  # chisq.scaling.factor.sb (added in 0.16); fall back to chisq.scaling.factor
  # with a warning for objects computed before it existed. The Satorra (2000)
  # method below does not use c0/c1 (it recomputes the traces), so this only
  # affects the SB 2001/2010 branches.
  sb_factor <- function(mm){
    sf <- mm@fitmeasures$chisq.scaling.factor.sb
    if (!is.null(sf) && length(sf) == 1 && is.finite(sf)) return(sf)
    sf2 <- mm@fitmeasures$chisq.scaling.factor
    if (!is.null(sf2) && length(sf2) == 1 && is.finite(sf2) &&
        !is.null(mm@fitmeasures$df) && mm@fitmeasures$df > 0){
      warning("This model was computed before the mean-adjusted scaling factor ",
              "(chisq.scaling.factor.sb) was stored; the Satorra-Bentler ",
              "difference test may be inaccurate for scaled-shifted / ",
              "mean-and-variance adjusted estimators. Re-run the model, or use ",
              "scaled.test.method = 'satorra.2000'.")
    }
    sf2
  }
  c1 <- sb_factor(m1)
  c0 <- sb_factor(m0)
  na_out <- list(stat = NA_real_, df = m, pvalue = NA_real_)

  if (!is.finite(m) || m <= 0) return(na_out)

  # Edge cases on the scaling factors (Satorra-Bentler 2001/2010) MUST be
  # applied before validating c0/c1: a saturated reference model (df 0) has no
  # scaling factor of its own (compute_ml_scaled_test returns NULL when df 0),
  # so its c is legitimately absent and is replaced by the convention below.
  #   r1 == 0 (M1 saturated): c1 is undefined -> use c1 = 1
  #   T1 ~ 0 with r1 > 0     : a perfect fit  -> use c1 = 0
  if (r1 == 0){
    c1 <- 1
  } else if (is.finite(T1) && abs(T1) < sqrt(.Machine$double.eps)){
    c1 <- 0
  }
  if (r0 == 0) c0 <- 1   # cannot happen for the larger-df model, but guard anyway

  if (is.null(c0) || is.null(c1) || !is.finite(c0) || !is.finite(c1) ||
      !is.finite(T0) || !is.finite(T1)) return(na_out)

  if (method == "satorra.2000"){
    A <- scaled_diff_A(m1, m0)
    if (is.null(A) || nrow(A) == 0){
      warning("Satorra (2000) scaled difference test could not be computed ",
              "(constraint Jacobian unavailable); reporting NA.")
      return(na_out)
    }
    tr <- scaled_diff_satorra2000_traces(m1, A)
    if (is.null(tr) || tr$trUG2 <= 0){
      warning("Satorra (2000) scaled difference test could not be computed ",
              "(trace assembly failed); reporting NA.")
      return(na_out)
    }
    if (scaled.shifted){
      a <- sqrt(m / tr$trUG2)
      b <- m - a * tr$trUG
      stat <- a * (T0 - T1) + b
    } else {
      cd <- tr$trUG / m
      if (!is.finite(cd) || cd <= 0){
        warning("Satorra (2000) scaled difference test has a non-positive ",
                "scaling factor; reporting NA.")
        return(na_out)
      }
      stat <- (T0 - T1) / cd
    }
  } else {
    if (method == "satorra.bentler.2010"){
      c10 <- scaled_diff_c10(m1, m0)
      if (is.null(c10) || !is.finite(c10)){
        warning("Satorra-Bentler (2010) scaled difference test could not be ",
                "computed (c10 unavailable); reporting NA.")
        return(na_out)
      }
      c1 <- c10
    }
    # Satorra-Bentler 2001 / 2010 mean-scaling:
    cd <- (r0 * c0 - r1 * c1) / m
    if (!is.finite(cd) || cd <= 0){
      warning("Scaled chi-square difference test has a non-positive scaling ",
              "factor (cd <= 0); reporting NA. Consider method = 'satorra.2000'.")
      return(na_out)
    }
    stat <- (T0 - T1) / cd
  }

  pval <- if (is.finite(stat)) pchisq(stat, m, lower.tail = FALSE) else NA_real_
  list(stat = stat, df = m, pvalue = pval)
}
