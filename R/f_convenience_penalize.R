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
    verbose,
    log = TRUE
) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")
  if (missing(verbose)) verbose <- x@verbose

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

    # Build the selection mask
    if (missing(row) && missing(col)) {
      # Penalize whole matrix (off-diagonal for symmetric, all for non-symmetric like beta)
      if (is_sym) {
        # Symmetric matrices: only off-diagonal (row != col)
        whichPen <- which(
          x@parameters$matrix == mat &
            !x@parameters$fixed &
            x@parameters$par > 0 &
            x@parameters$row != x@parameters$col &
            x@parameters$group_id %in% group
        )
      } else {
        # Non-symmetric matrices (e.g., beta): all elements
        whichPen <- which(
          x@parameters$matrix == mat &
            !x@parameters$fixed &
            x@parameters$par > 0 &
            x@parameters$group_id %in% group
        )
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

      whichPen <- which(
        x@parameters$matrix == mat &
          x@parameters$row %in% r &
          x@parameters$col %in% cl &
          !x@parameters$fixed &
          x@parameters$par > 0 &
          x@parameters$group_id %in% group
      )
    }

    # Set penalty_lambda for matching rows:
    if (length(whichPen) > 0) {
      x@parameters$penalty_lambda[whichPen] <- lambda
      nPenalized <- nPenalized + length(whichPen)
    }
  }

  if (verbose) {
    message(paste0("Penalized ", nPenalized, " parameter table entries (lambda = ", lambda, ")"))
  }

  if (log) {
    x <- addLog(x, paste0("Penalized ", nPenalized, " entries in matrix/matrices: ",
                           paste(matrix, collapse = ", "), " (lambda = ", lambda, ")"))
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
