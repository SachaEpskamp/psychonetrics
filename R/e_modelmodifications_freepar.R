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
  startEPC = TRUE,
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
  } else {
    expected <-  x@parameters$est[whichFree][!is.na(x@parameters$epc[whichFree])]      +  x@parameters$epc[whichFree][!is.na(x@parameters$epc[whichFree])]    
    if (startEPC){
      # Set to EPC:
      x@parameters$est[whichFree][!is.na(x@parameters$epc[whichFree])] <-   expected
    } else {
      # Set to EPC:
      x@parameters$est[whichFree][!is.na(x@parameters$epc[whichFree])] <- 0.0001*sign(expected)
    }
  }
  
  # Adjust start values in the whole matrix?
  if (!startEPC){
    inds <- x@parameters$matrix == x@parameters$matrix[whichFree] & !x@parameters$fixed & !x@parameters$identified
    x@parameters$est[inds] <- 
      0.0001*sign(x@parameters$est[inds])
  }

  
  # x@parameters$std[whichFree] <- NA
  x@parameters$par[whichFree] <- curMax + seq_len(length(whichFree))
  # x@parameters$se[whichFree] <- NA
  # x@parameters$p[whichFree] <- NA
  # x@parameters$mi[whichFree] <- NA
  # x@parameters$pmi[whichFree] <- NA
  # x@parameters$mi_equal[whichFree] <- NA
  # x@parameters$pmi_equal[whichFree] <- NA
  x@parameters$fixed[whichFree] <- FALSE
  
  x@parameters <- clearpars(x@parameters,whichFree)

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