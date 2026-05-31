# Function to constrain all (free) elements of a matrix to be equal WITHIN each
# group. In contrast to groupequal() -- which constrains corresponding elements
# equal ACROSS groups -- allequal() shares a single parameter across all selected
# elements of the matrix within a group, and does so separately for each group
# (so the constraint is not imposed across groups). This is convenient for, e.g.,
# constraining all thresholds (tau) or all Blume-Capel quadratic potentials
# (delta) to a common value, or fitting an equal-edge-weight network (omega).
allequal <- function(
  x, # Model
  matrix, # Must be given
  row, # Can be missing (default: all rows of the matrix)
  col, # Can be missing (default: all columns of the matrix)
  group, # Can be missing (default: all groups)
  verbose,
  log = TRUE,
  runmodel = FALSE,
  identify = TRUE,
  ...){

  if (missing(verbose)){
    verbose <- x@verbose
  }

  if (missing(matrix)){
    stop("'matrix' argument may not be missing")
  }
  if (length(matrix) > 1){
    stop("'matrix' must be a string of length 1")
  }
  if (!matrix %in% x@parameters$matrix){
    stop(paste0("'matrix' argument must be one of the modeled matrices: ", paste0(unique(x@matrices$name),collapse=", ")))
  }

  # Default to all groups:
  if (missing(group)){
    group <- x@sample@groups$id
  }
  # Default to all rows/cols of the matrix:
  if (missing(row)){
    row <- unique(x@parameters$row[x@parameters$matrix == matrix])
  }
  if (missing(col)){
    col <- unique(x@parameters$col[x@parameters$matrix == matrix])
  }
  # If row/col is character, convert to number:
  if (is.character(row) | is.character(col)){
    labs <- labtoind(x,row,col,matrix)
  }
  if (is.character(row)){
    row <- labs$row
  }
  if (is.character(col)){
    col <- labs$col
  }
  # If the matrix is symmetric, make sure col contains the lower number:
  if (x@matrices$symmetric[x@matrices$name == matrix]){
    row0 <- row
    col0 <- col

    col <- pmin(row0,col0)
    row <- pmax(row0,col0)
  }

  # current max par:
  curMax <- max(x@parameters$par)

  # Set the model to be not computed:
  x@computed <- FALSE

  # For each group, set all selected free elements equal (share one parameter).
  # The shared parameter differs between groups, so equality is not imposed
  # across groups (use groupequal() in addition for that):
  for (g in group){
    inds <- which(
      x@parameters$matrix == matrix &
      x@parameters$row %in% row &
      x@parameters$col %in% col &
      !x@parameters$fixed &
      x@parameters$group_id == g)

    if (length(inds) > 1){
      x@parameters$par[inds] <- x@parameters$par[inds[1]]
      # Use the mean of the current values as the common starting value:
      x@parameters$est[inds] <- mean(x@parameters$est[inds], na.rm = TRUE)
      x@parameters <- clearpars(x@parameters, inds)
    }
  }

  # Relabel:
  x@parameters <- parRelabel(x@parameters)

  # Identify:
  if (identify){
    x <- identify(x)
  }

  # Number of parameters constrained:
  nConstrained <- curMax - max(x@parameters$par)

  # Output:
  if (verbose){
    if (nConstrained == 0){
      message("No parameters needed to be constrained...")
    } else {
      message(paste0("Constrained ", nConstrained, " parameters!"))
    }
  }

  # Write to log:
  if (log){
    x <- addLog(x, paste0("Set all element(s) of ", matrix, " equal within groups. Constrained ", nConstrained, " parameters!"))
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }

  return(x)
}
