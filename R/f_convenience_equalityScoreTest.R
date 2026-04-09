## equalityScoreTest: score (Lagrange multiplier) test for cross-group equality
## constraints in a psychonetrics model. Inspired by lavaan::lavTestScore.
##
## Returns (invisibly) a list with two data frames:
##   $uni   - per released parameter (one row per (matrix,row,col,group)),
##            sorted by X2 descending. Mirrors the "univariate" output of
##            lavaan::lavTestScore.
##   $total - per equality-constrained (matrix,row,col), the joint
##            multivariate score test (df = G-1 typically), sorted by X2.

# ---------------------------------------------------------------------------
# Inner helper. Builds an augmented model in which all currently equality-
# constrained parameters are released across groups (G-1 new parameters per
# constrained tuple), computes the gradient and Fisher information, and from
# them the score statistics. Avoids the additional "free all zero parameters"
# augmentation done by addMIs_inner_full(type = "free"), so the resulting
# joint statistics are uncontaminated.
# ---------------------------------------------------------------------------
.equalityScoreTestInner <- function(x, matrices = NULL, analyticFisher = TRUE,
                                    method = c("jacobian","schur")){
  method <- match.arg(method)
  if (is.null(matrices)) matrices <- x@matrices$name

  parTab <- x@parameters
  G <- nrow(x@sample@groups)
  if (G < 2) return(list(uni = NULL, total = NULL))

  # Identify equality-constrained (matrix,row,col) tuples in the requested matrices.
  # An equality-constrained tuple has G rows, all with the SAME non-zero par.
  eq_set <- parTab %>%
    filter(.data[["par"]] != 0, .data[["matrix"]] %in% matrices,
           !.data[["identified"]], !.data[["fixed"]]) %>%
    group_by(.data[["matrix"]], .data[["row"]], .data[["col"]]) %>%
    summarize(npar = dplyr::n_distinct(.data[["par"]]),
              ngrp = dplyr::n(),
              .groups = "drop") %>%
    filter(.data[["npar"]] == 1, .data[["ngrp"]] == G)

  if (nrow(eq_set) == 0){
    return(list(uni = NULL, total = NULL))
  }

  # Build augmented copy: only release equality constraints, nothing else.
  modCopy <- x
  curMax <- max(modCopy@parameters$par)
  next_par <- curMax

  # Per-tuple metadata to map new pars back to (matrix,row,col,group):
  tuples <- vector("list", nrow(eq_set))
  for (i in seq_len(nrow(eq_set))){
    mat_i <- eq_set$matrix[i]; row_i <- eq_set$row[i]; col_i <- eq_set$col[i]
    rows_idx <- which(modCopy@parameters$matrix == mat_i &
                        modCopy@parameters$row == row_i &
                        modCopy@parameters$col == col_i)
    rows_idx <- rows_idx[order(modCopy@parameters$group_id[rows_idx])]
    k_cur <- modCopy@parameters$par[rows_idx[1]]  # group-1 (reference) par index
    new_pars <- integer(0)
    for (j in 2:length(rows_idx)){
      next_par <- next_par + 1L
      modCopy@parameters$par[rows_idx[j]] <- next_par
      new_pars <- c(new_pars, next_par)
    }
    tuples[[i]] <- list(
      mat = mat_i, row = row_i, col = col_i,
      var1 = modCopy@parameters$var1[rows_idx[1]],
      op   = modCopy@parameters$op[rows_idx[1]],
      var2 = modCopy@parameters$var2[rows_idx[1]],
      k_cur = k_cur,
      new_pars = new_pars,
      group_ids = modCopy@parameters$group_id[rows_idx[-1]],
      group_labels = modCopy@parameters$group[rows_idx[-1]]
    )
  }

  # Rebuild M matrix:
  if (modCopy@cpp){
    modCopy@extramatrices$M <- Mmatrix_cpp(modCopy@parameters)
  } else {
    modCopy@extramatrices$M <- Mmatrix(modCopy@parameters)
  }

  # Gradient and Fisher information at the (restricted) MLE:
  if (modCopy@cpp){
    g <- psychonetrics_gradient_cpp(parVector(modCopy), modCopy)
  } else {
    g <- psychonetrics_gradient(parVector(modCopy), modCopy)
  }
  nTotal <- sum(x@sample@groups$nobs)
  if (modCopy@cpp){
    H <- 4 * nTotal * as(psychonetrics_FisherInformation_cpp(modCopy, analyticFisher), "matrix")
  } else {
    H <- 4 * nTotal * as(psychonetrics_FisherInformation(modCopy, analyticFisher), "matrix")
  }

  curInds <- seq_len(curMax)
  newInds <- (curMax + 1L):next_par
  if (length(newInds) == 0L){
    return(list(uni = NULL, total = NULL))
  }
  uni_rows <- vector("list", 0)
  total_rows <- vector("list", 0)

  if (method == "jacobian"){
    # Reproduce lavaan::lavTestScore exactly. Translation of psychonetrics
    # quantities to lavaan's:
    #   score_lav = -(N/2) * g       (lavaan's gradient.logl = dlogL/dθ)
    #   Info_lav  = H / (8 * N)      (lavaan's information.expected, per obs)
    # because psychonetrics works in F = -(2/N) logL, so dlogL/dθ = -(N/2) g
    # and J_per_obs = I_F / 2 = (H/(4N))/2 = H/(8N).
    npar <- ncol(H)
    score_lav <- -(nTotal / 2) * g
    Info_lav  <- H / (8 * nTotal)

    # Build the constraint Jacobian R: one row per released parameter.
    # Row r has +1 at the reference (k_cur) column and -1 at the released
    # (k_new) column. We also remember which tuple each row belongs to so
    # the per-tuple joint test can release G-1 rows simultaneously.
    R_list <- list()
    row_tuple <- integer(0)
    row_within_tuple <- integer(0)
    for (i in seq_along(tuples)){
      info <- tuples[[i]]
      for (j in seq_along(info$new_pars)){
        rr <- numeric(npar)
        rr[info$k_cur] <- 1
        rr[info$new_pars[j]] <- -1
        R_list[[length(R_list) + 1L]] <- rr
        row_tuple <- c(row_tuple, i)
        row_within_tuple <- c(row_within_tuple, j)
      }
    }
    R <- do.call(rbind, R_list)

    # Helper closure: compute Z1.plus1 = upper-left npar block of
    # ginv([Info_lav, R1'; R1, 0]) for a given R1.
    Z1plus1 <- function(R1){
      if (nrow(R1) == 0L){
        return(MASS::ginv(Info_lav))
      }
      Z1 <- rbind(
        cbind(Info_lav, t(R1)),
        cbind(R1, matrix(0, nrow(R1), nrow(R1)))
      )
      Z1p <- MASS::ginv(Z1)
      Z1p[seq_len(npar), seq_len(npar)]
    }

    # Univariate: one stat per row of R (release just that one row).
    for (r in seq_len(nrow(R))){
      R1 <- R[-r, , drop = FALSE]
      Z <- tryCatch(Z1plus1(R1), error = function(e) NULL)
      X2_uni <- if (is.null(Z)) NA_real_
      else as.numeric(t(score_lav) %*% Z %*% score_lav / (2 * nTotal))
      tup_i <- row_tuple[r]; j_in_tup <- row_within_tuple[r]
      info <- tuples[[tup_i]]
      uni_rows[[length(uni_rows) + 1L]] <- data.frame(
        var1 = info$var1, op = info$op, var2 = info$var2,
        matrix = info$mat, row = info$row, col = info$col,
        group = info$group_labels[j_in_tup], group_id = info$group_ids[j_in_tup],
        X2 = X2_uni, df = 1L,
        p.value = if (is.na(X2_uni)) NA_real_ else pchisq(X2_uni, 1, lower.tail = FALSE),
        stringsAsFactors = FALSE
      )
    }

    # Joint per tuple: release the G-1 rows belonging to that tuple at once.
    for (i in seq_along(tuples)){
      info <- tuples[[i]]
      r_set <- which(row_tuple == i)
      R1 <- R[-r_set, , drop = FALSE]
      Z <- tryCatch(Z1plus1(R1), error = function(e) NULL)
      Tjoint <- if (is.null(Z)) NA_real_
      else as.numeric(t(score_lav) %*% Z %*% score_lav / (2 * nTotal))
      df_j <- length(r_set)
      p_j <- if (is.na(Tjoint)) NA_real_ else pchisq(Tjoint, df_j, lower.tail = FALSE)
      total_rows[[length(total_rows) + 1L]] <- data.frame(
        var1 = info$var1, op = info$op, var2 = info$var2,
        matrix = info$mat, row = info$row, col = info$col,
        X2 = Tjoint, df = df_j, p.value = p_j,
        stringsAsFactors = FALSE
      )
    }
  } else {
    # Schur-complement fallback (faster for complex models). NOT exactly
    # equivalent to lavaan::lavTestScore in unbalanced multi-group samples,
    # but matches in balanced cases. Computes T = 2 * (n*g_S)' V[S,S]^{-1} (n*g_S)
    # with V the Schur complement of the augmented Fisher information.
    V <- H[newInds, newInds, drop = FALSE] -
      H[newInds, curInds, drop = FALSE] %*%
      solve_symmetric(H[curInds, curInds]) %*%
      H[curInds, newInds, drop = FALSE]
    for (i in seq_along(tuples)){
      info <- tuples[[i]]
      S_par <- info$new_pars
      S_pos <- S_par - curMax
      Vsub <- V[S_pos, S_pos, drop = FALSE]
      gsub <- (-nTotal) * g[S_par]
      Vdiag <- diag(Vsub)
      for (j in seq_along(S_pos)){
        X2_uni <- if (!is.finite(Vdiag[j]) || Vdiag[j] < sqrt(.Machine$double.eps)) {
          NA_real_
        } else {
          2 * (gsub[j]^2) / Vdiag[j]
        }
        uni_rows[[length(uni_rows) + 1L]] <- data.frame(
          var1 = info$var1, op = info$op, var2 = info$var2,
          matrix = info$mat, row = info$row, col = info$col,
          group = info$group_labels[j], group_id = info$group_ids[j],
          X2 = X2_uni, df = 1L,
          p.value = if (is.na(X2_uni)) NA_real_ else pchisq(X2_uni, 1, lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }
      Tjoint <- tryCatch(
        as.numeric(2 * t(gsub) %*% solve(Vsub) %*% gsub),
        error = function(e) NA_real_
      )
      df_j <- length(S_pos)
      p_j <- if (is.na(Tjoint)) NA_real_ else pchisq(Tjoint, df_j, lower.tail = FALSE)
      total_rows[[length(total_rows) + 1L]] <- data.frame(
        var1 = info$var1, op = info$op, var2 = info$var2,
        matrix = info$mat, row = info$row, col = info$col,
        X2 = Tjoint, df = df_j, p.value = p_j,
        stringsAsFactors = FALSE
      )
    }
  }

  uni <- if (length(uni_rows) > 0) do.call(rbind, uni_rows) else NULL
  total <- if (length(total_rows) > 0) do.call(rbind, total_rows) else NULL
  if (!is.null(uni)){
    uni <- uni[order(uni$X2, decreasing = TRUE),]
    rownames(uni) <- NULL
  }
  if (!is.null(total)){
    total <- total[order(total$X2, decreasing = TRUE),]
    rownames(total) <- NULL
  }
  list(uni = uni, total = total)
}


# ---------------------------------------------------------------------------
# User-facing function.
# ---------------------------------------------------------------------------
equalityScoreTest <- function(x, matrices, top = 10, verbose = TRUE,
                              method = c("jacobian","schur")){
  method <- match.arg(method)
  if (!is(x, "psychonetrics")){
    stop("'x' must be a psychonetrics model.")
  }
  if (nrow(x@sample@groups) < 2){
    stop("equalityScoreTest() requires a multi-group model.")
  }
  if (!x@computed){
    stop("Model has not been run. Use runmodel() first.")
  }
  if (missing(matrices)) matrices <- x@matrices$name

  res <- .equalityScoreTestInner(x, matrices = matrices, method = method)

  if (verbose){
    cat("\nUnivariate score tests (one per released parameter):\n\n")
    if (is.null(res$uni) || nrow(res$uni) == 0){
      cat("  (no equality-constrained parameters found)\n")
    } else {
      uni_print <- res$uni
      uni_print$X2 <- goodNum(uni_print$X2)
      uni_print$p.value <- goodNum(uni_print$p.value)
      print.data.frame(head(uni_print, top), row.names = FALSE)
    }

    cat("\nJoint score tests (one per equality-constrained parameter, df = G-1):\n\n")
    if (is.null(res$total) || nrow(res$total) == 0){
      cat("  (no equality-constrained parameters found)\n")
    } else {
      total_print <- res$total
      total_print$X2 <- goodNum(total_print$X2)
      total_print$p.value <- goodNum(total_print$p.value)
      print.data.frame(head(total_print, top), row.names = FALSE)
    }
  }

  invisible(res)
}
