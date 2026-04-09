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
.equalityScoreTestInner <- function(x, matrices = NULL, analyticFisher = TRUE){
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
  V <- H[newInds, newInds, drop = FALSE] -
    H[newInds, curInds, drop = FALSE] %*%
    solve_symmetric(H[curInds, curInds]) %*%
    H[curInds, newInds, drop = FALSE]

  uni_rows <- vector("list", 0)
  total_rows <- vector("list", 0)
  for (i in seq_along(tuples)){
    info <- tuples[[i]]
    S_par <- info$new_pars
    S_pos <- S_par - curMax  # position within newInds
    Vsub <- V[S_pos, S_pos, drop = FALSE]
    gsub <- (-nTotal) * g[S_par]
    # Univariate per released parameter (diagonal of V).
    # The factor of 2 below converts the psychonetrics discrepancy-function-based
    # statistic to the -2 log L scale used by lavaan::lavTestScore, so that the
    # values are directly comparable to lavaan's output.
    Vdiag <- diag(Vsub)
    for (k in seq_along(S_pos)){
      X2_uni <- if (!is.finite(Vdiag[k]) || Vdiag[k] < sqrt(.Machine$double.eps)) {
        NA_real_
      } else {
        2 * (gsub[k]^2) / Vdiag[k]
      }
      uni_rows[[length(uni_rows) + 1L]] <- data.frame(
        var1 = info$var1, op = info$op, var2 = info$var2,
        matrix = info$mat, row = info$row, col = info$col,
        group = info$group_labels[k], group_id = info$group_ids[k],
        X2 = X2_uni, df = 1L,
        p.value = if (is.na(X2_uni)) NA_real_ else pchisq(X2_uni, 1, lower.tail = FALSE),
        stringsAsFactors = FALSE
      )
    }
    # Joint stat (factor 2 to match lavaan -2 log L scale):
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
equalityScoreTest <- function(x, matrices, top = 10, verbose = TRUE){
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

  res <- .equalityScoreTestInner(x, matrices = matrices)

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
