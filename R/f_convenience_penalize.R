# Functions to penalize / unpenalize parameters for PML estimation
# Follows the pattern of fixpar() / freepar()

# Helper to extract per-free-parameter penalty lambda from the parameter table
penaltyVector <- function(x) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")
  nPar <- max(x@parameters$par)
  if (nPar == 0) return(numeric(0))
  lambdas <- numeric(nPar)
  for (i in seq_len(nPar)) {
    rows <- which(x@parameters$par == i)
    if (length(rows) > 0) {
      val <- x@parameters$penalty_lambda[rows[1]]
      # NA means auto-select (not yet resolved); treat as 0 for optimizer safety
      lambdas[i] <- if (is.na(val)) 0 else val
    }
  }
  lambdas
}

# Check if model has any auto-select (NA) penalty_lambda values
hasAutoLambda <- function(x) {
  any(is.na(x@parameters$penalty_lambda) & x@parameters$par > 0 & !x@parameters$fixed)
}


# Penalize: set penalty_lambda for matching parameters
penalize <- function(
    x,            # psychonetrics model
    matrix,       # Character: matrix name(s) to penalize. Default: defaultPenalizeMatrices(x)
    row,          # Optional: specific rows (if missing, whole matrix)
    col,          # Optional: specific cols
    lambda,       # Optional: per-parameter lambda. Default: x@penalty$lambda
    group,        # Optional: specific group(s)
    what = c("free", "fixed"),  # Penalize free or fixed parameters
    verbose,
    log = TRUE
) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")
  if (missing(verbose)) verbose <- x@verbose
  what <- match.arg(what)

  # Default matrix selection:
  if (missing(matrix)) {
    matrix <- defaultPenalizeMatrices(x)
  }

  # Default lambda from global penalty:
  if (missing(lambda)) {
    lambda <- x@penalty$lambda
    if (is.null(lambda) || (!is.na(lambda) && lambda == 0)) {
      warning("No lambda specified and global penalty lambda is 0. Set lambda > 0 for penalization.")
    }
  }

  # Default to all groups:
  if (missing(group)) {
    group <- x@sample@groups$id
  }

  # Track how many parameters were penalized:
  nPenalized <- 0

  # Check if lambda is a matrix (per-element penalty specification):
  lambda_is_matrix <- is.matrix(lambda) || (is.array(lambda) && length(dim(lambda)) >= 2)

  for (mat in matrix) {
    # Validate matrix name:
    if (!mat %in% x@parameters$matrix) {
      warning(paste0("Matrix '", mat, "' not found in model. Skipping."))
      next
    }

    # Is this matrix symmetric?
    mat_info <- x@matrices[x@matrices$name == mat, ]
    is_sym <- FALSE
    if (nrow(mat_info) > 0) {
      is_sym <- mat_info$symmetrical[1]
    }

    if (lambda_is_matrix) {
      # --- Matrix-form lambda: per-element penalty specification ---
      lambda_mat <- if (length(dim(lambda)) == 3) lambda[,,1] else as.matrix(lambda)

      # Validate dimensions:
      mat_nrow <- max(x@parameters$row[x@parameters$matrix == mat])
      mat_ncol <- max(x@parameters$col[x@parameters$matrix == mat])
      if (nrow(lambda_mat) != mat_nrow || ncol(lambda_mat) != mat_ncol) {
        stop("Dimensions of lambda matrix (", nrow(lambda_mat), "x", ncol(lambda_mat),
             ") don't match model matrix '", mat, "' (", mat_nrow, "x", mat_ncol, ")")
      }

      nFreed <- 0
      for (i in seq_len(nrow(lambda_mat))) {
        for (j in seq_len(ncol(lambda_mat))) {
          # For symmetric matrices, only process lower triangle (including diagonal):
          if (is_sym && i < j) next

          val <- lambda_mat[i, j]

          # Determine row/col in parameter table (symmetric: row >= col):
          r <- i; cl <- j
          if (is_sym) { r <- max(i, j); cl <- min(i, j) }

          # Find matching parameter table rows:
          matching <- which(
            x@parameters$matrix == mat &
              x@parameters$row == r &
              x@parameters$col == cl &
              x@parameters$group_id %in% group
          )
          if (length(matching) == 0) next

          if (!is.na(val) && val == 0) {
            # 0 = don't penalize (set penalty_lambda to 0)
            x@parameters$penalty_lambda[matching] <- 0
          } else {
            # NA or > 0 = penalize this element
            if (what == "free") {
              # Only penalize already-free parameters (silently skip fixed):
              matching_free <- matching[!x@parameters$fixed[matching] &
                                         x@parameters$par[matching] > 0]
              if (length(matching_free) > 0) {
                x@parameters$penalty_lambda[matching_free] <- val
                nPenalized <- nPenalized + length(matching_free)
              }
            } else {
              # what == "fixed": free fixed parameters and penalize them (silently skip free):
              for (idx in matching) {
                if (x@parameters$fixed[idx] && x@parameters$par[idx] == 0) {
                  curMax <- max(x@parameters$par)
                  x@parameters$par[idx] <- curMax + 1
                  x@parameters$fixed[idx] <- FALSE
                  x@parameters$est[idx] <- 0  # Start at 0 for penalized params
                  nFreed <- nFreed + 1
                  # Set penalty:
                  x@parameters$penalty_lambda[idx] <- val
                  nPenalized <- nPenalized + 1
                }
              }
            }
          }
        }
      }

      # Relabel parameter indices after freeing:
      if (nFreed > 0) {
        x@parameters <- parRelabel(x@parameters)
      }

    } else {
      # --- Scalar lambda: existing behavior ---

      if (missing(row) && missing(col)) {
        # Penalize whole matrix
        if (what == "free") {
          # Penalize free parameters (current default behavior)
          if (is_sym) {
            whichPen <- which(
              x@parameters$matrix == mat &
                !x@parameters$fixed &
                x@parameters$par > 0 &
                x@parameters$row != x@parameters$col &
                x@parameters$group_id %in% group
            )
          } else {
            whichPen <- which(
              x@parameters$matrix == mat &
                !x@parameters$fixed &
                x@parameters$par > 0 &
                x@parameters$group_id %in% group
            )
          }
        } else {
          # what == "fixed": free fixed parameters and penalize them
          if (is_sym) {
            whichFixed <- which(
              x@parameters$matrix == mat &
                x@parameters$fixed &
                x@parameters$par == 0 &
                x@parameters$row != x@parameters$col &
                x@parameters$group_id %in% group
            )
          } else {
            whichFixed <- which(
              x@parameters$matrix == mat &
                x@parameters$fixed &
                x@parameters$par == 0 &
                x@parameters$group_id %in% group
            )
          }
          # Free them:
          if (length(whichFixed) > 0) {
            for (idx in whichFixed) {
              curMax <- max(x@parameters$par)
              x@parameters$par[idx] <- curMax + 1
              x@parameters$fixed[idx] <- FALSE
              x@parameters$est[idx] <- 0  # Start at 0 for penalized params
            }
            x@parameters <- parRelabel(x@parameters)
          }
          whichPen <- whichFixed
        }
      } else {
        # Specific row/col selection
        r <- if (!missing(row)) row else unique(x@parameters$row[x@parameters$matrix == mat])
        cl <- if (!missing(col)) col else unique(x@parameters$col[x@parameters$matrix == mat])

        # If row/col are character, convert to indices
        if (is.character(r) || is.character(cl)) {
          labs <- labtoind(x, r, cl, mat)
          if (is.character(r)) r <- labs$row
          if (is.character(cl)) cl <- labs$col
        }

        # For symmetric matrices, ensure consistent ordering (row >= col)
        if (is_sym) {
          r0 <- r; cl0 <- cl
          r <- pmax(r0, cl0)
          cl <- pmin(r0, cl0)
        }

        # Find all matching entries:
        allMatching <- which(
          x@parameters$matrix == mat &
            x@parameters$row %in% r &
            x@parameters$col %in% cl &
            x@parameters$group_id %in% group
        )

        if (what == "free") {
          # Penalize free parameters; warn about fixed ones:
          whichPen <- allMatching[!x@parameters$fixed[allMatching] &
                                    x@parameters$par[allMatching] > 0]
          whichSkipped <- allMatching[x@parameters$fixed[allMatching] |
                                       x@parameters$par[allMatching] == 0]
          if (length(whichSkipped) > 0) {
            warning("Cannot penalize ", length(whichSkipped), " fixed parameter(s) in '", mat,
                    "' when what = 'free'. Use what = 'fixed' to first free these parameters before penalizing them.")
          }
        } else {
          # what == "fixed": free fixed parameters and penalize; warn about free ones:
          whichFixed <- allMatching[x@parameters$fixed[allMatching] &
                                     x@parameters$par[allMatching] == 0]
          whichSkipped <- allMatching[!x@parameters$fixed[allMatching] &
                                       x@parameters$par[allMatching] > 0]
          if (length(whichSkipped) > 0) {
            warning("Cannot penalize ", length(whichSkipped), " already-free parameter(s) in '", mat,
                    "' when what = 'fixed'. Use what = 'free' to penalize parameters that are already free.")
          }
          # Free them:
          if (length(whichFixed) > 0) {
            for (idx in whichFixed) {
              curMax <- max(x@parameters$par)
              x@parameters$par[idx] <- curMax + 1
              x@parameters$fixed[idx] <- FALSE
              x@parameters$est[idx] <- 0  # Start at 0 for penalized params
            }
            x@parameters <- parRelabel(x@parameters)
          }
          whichPen <- whichFixed
        }
      }

      # Set penalty_lambda for matching rows:
      if (length(whichPen) > 0) {
        x@parameters$penalty_lambda[whichPen] <- lambda
        nPenalized <- nPenalized + length(whichPen)
      }
    }
  }

  if (verbose) {
    if (lambda_is_matrix) {
      message(paste0("Penalized ", nPenalized, " parameter table entries (matrix-form lambda, what = '", what, "')"))
    } else {
      message(paste0("Penalized ", nPenalized, " parameter table entries (lambda = ", lambda, ", what = '", what, "')"))
    }
  }

  if (log) {
    lambda_str <- if (lambda_is_matrix) "matrix" else as.character(lambda)
    x <- addLog(x, paste0("Penalized ", nPenalized, " entries in matrix/matrices: ",
                           paste(matrix, collapse = ", "), " (lambda = ", lambda_str, ", what = '", what, "')"))
  }

  # Set model to not computed:
  x@computed <- FALSE

  x
}


# Unpenalize: remove penalty from matching parameters
unpenalize <- function(
    x,            # psychonetrics model
    matrix,       # Character: matrix name(s) to unpenalize
    row,          # Optional: specific rows
    col,          # Optional: specific cols
    group,        # Optional: specific group(s)
    verbose,
    log = TRUE
) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")
  if (missing(verbose)) verbose <- x@verbose

  # Default matrix selection:
  if (missing(matrix)) {
    matrix <- unique(x@parameters$matrix[!is.na(x@parameters$penalty_lambda) & x@parameters$penalty_lambda > 0])
  }

  # Default to all groups:
  if (missing(group)) {
    group <- x@sample@groups$id
  }

  nUnpenalized <- 0

  for (mat in matrix) {
    if (!mat %in% x@parameters$matrix) {
      warning(paste0("Matrix '", mat, "' not found in model. Skipping."))
      next
    }

    mat_info <- x@matrices[x@matrices$name == mat, ]
    is_sym <- FALSE
    if (nrow(mat_info) > 0) is_sym <- mat_info$symmetrical[1]

    if (missing(row) && missing(col)) {
      whichUnpen <- which(
        x@parameters$matrix == mat &
          !is.na(x@parameters$penalty_lambda) &
          x@parameters$penalty_lambda > 0 &
          x@parameters$group_id %in% group
      )
    } else {
      r <- if (!missing(row)) row else unique(x@parameters$row[x@parameters$matrix == mat])
      cl <- if (!missing(col)) col else unique(x@parameters$col[x@parameters$matrix == mat])

      if (is.character(r) || is.character(cl)) {
        labs <- labtoind(x, r, cl, mat)
        if (is.character(r)) r <- labs$row
        if (is.character(cl)) cl <- labs$col
      }

      if (is_sym) {
        r0 <- r; cl0 <- cl
        r <- pmax(r0, cl0)
        cl <- pmin(r0, cl0)
      }

      whichUnpen <- which(
        x@parameters$matrix == mat &
          x@parameters$row %in% r &
          x@parameters$col %in% cl &
          !is.na(x@parameters$penalty_lambda) &
          x@parameters$penalty_lambda > 0 &
          x@parameters$group_id %in% group
      )
    }

    if (length(whichUnpen) > 0) {
      x@parameters$penalty_lambda[whichUnpen] <- 0
      nUnpenalized <- nUnpenalized + length(whichUnpen)
    }
  }

  if (verbose) {
    message(paste0("Unpenalized ", nUnpenalized, " parameter table entries"))
  }

  if (log) {
    x <- addLog(x, paste0("Unpenalized ", nUnpenalized, " entries in matrix/matrices: ",
                           paste(matrix, collapse = ", ")))
  }

  x@computed <- FALSE
  x
}
