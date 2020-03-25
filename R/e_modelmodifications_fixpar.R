# Function to fix a certain parameter
fixpar <- function(
  x, # Model
  matrix, # Must be given
  row, # Must be given
  col, # Must be given
  value = 0, # typical so can be missing
  group, # Can be missing to indicate all
  verbose,
  log = TRUE,
  runmodel = FALSE,
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
  # If row is character, convert to number:
  if (is.character(row) | is.character(col)){
    labs <- labtoind(x,row,col,matrix)
  }
  if (is.character(row)){
    row <- labs$row
  }
  if (is.character(col)){
    col <- labs$col
  }
  # If the matrix is symmetric, add them to each other:
  if (x@matrices$symmetric[x@matrices$name == matrix]){
    row0 <- row
    col0 <- col
    row <- c(row0,col0)
    col <- c(col0,row0)
  }
  
  # which to fix:
  whichCons <- which(x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & !x@parameters$fixed & x@parameters$group_id %in% group)
  
  # Length0?
  if (length(whichCons) == 0){
    if (verbose) message("No parameters need to be fixed")
      return(x)
  }
  
  # current max par:
  curMax <- max(x@parameters$par)
  
  # Set the model to be not computed:
  x@computed <- FALSE
  
  # Fix the parameters:
  x@parameters$est[whichCons] <- value
  # x@parameters$std[whichCons] <- NA
  x@parameters$par[whichCons] <- 0
  # x@parameters$se[whichCons] <- NA
  # x@parameters$p[whichCons] <- NA
  # x@parameters$mi[whichCons] <- NA
  # x@parameters$pmi[whichCons] <- NA
  # x@parameters$mi_equal[whichCons] <- NA
  # x@parameters$pmi_equal[whichCons] <- NA
  x@parameters$fixed[whichCons] <- TRUE
  x@parameters <- clearpars(x@parameters,whichCons)

  # Relabel:
  x@parameters   <- parRelabel(x@parameters)
  
  # Output:
  if (verbose){
    message(paste0("Fixed ",curMax - max(x@parameters$par)," parameters!"))
  }
  

  # Write to log:
  if (log){
    # Add log:
    x <- addLog(x, paste0("Fixed element(s) of ",matrix,": ",curMax - max(x@parameters$par)," parameters!")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  

  return(x)
}