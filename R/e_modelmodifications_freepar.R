# Function to fix a certain parameter
freepar <- function(
  x, # Model
  matrix, # Must be given
  row, # Must be given
  col, # Must be given
  start, # Starting value, can be ignored
  group, # Can be missing to indicate all
  verbose = TRUE,
  log = TRUE,
  runmodel = FALSE,
  ...){

  if (missing(matrix)){
    stop("'matrix' argument may not be missing")
  }
  if (length(matrix) > 1){
    stop("'matrix' must be a string of length 1")
  }
  if (!matrix %in% x@parameters$matrix){
    stop(paste0("'matrix' argument must be one of the modeled matrices: ", paste0(unique(x@matrices$name),collapse=", ")))
  }
  if (missing(row)){
    stop("'row' argument must be given")
  }
  if (missing(col)){
    stop("'col' argument must be given")
  }
  
  # If groups is missing, just do all groups:
  if (missing(group)){
    group <- x@sample@groups$id
  }

  # If the matrix is symmetric, add them to each other:
  if (x@matrices$symmetric[x@matrices$name == matrix]){
    row0 <- row
    col0 <- col
    row <- c(row0,col0)
    col <- c(col0,row0)
  }

  # which to free:
  whichFree <- which(x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & x@parameters$fixed & x@parameters$group_id %in% group)
  
  # Length0?
  if (length(whichFree) == 0){
    if (verbose) message("No parameters need to be freed")
      return(x)
  }
  
  # current max par:
  curMax <- max(x@parameters$par)
  
  # Set the model to be not computed:
  x@computed <- FALSE
  
  # Fix the parameters:
  if (!missing(start)){
    x@parameters$est[whichFree] <- start  
  }
  x@parameters$std[whichFree] <- NA
  x@parameters$par[whichFree] <- curMax + seq_len(length(whichFree))
  x@parameters$se[whichFree] <- NA
  x@parameters$p[whichFree] <- NA
  x@parameters$mi[whichFree] <- NA
  x@parameters$pmi[whichFree] <- NA
  x@parameters$mi_equal[whichFree] <- NA
  x@parameters$pmi_equal[whichFree] <- NA
  x@parameters$fixed[whichFree] <- FALSE

  # Relabel:
  x@parameters   <- parRelabel(x@parameters)
  
  # Output:
  if (verbose){
    message(paste0("Freed ",max(x@parameters$par) - curMax," parameters!"))
  }
  

  # Write to log:
  if (log){
    # Add log:
    x <- addLog(x, paste0("Freed element(s) of ",matrix,": ",max(x@parameters$par) - curMax," parameters!")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  

  return(x)
}