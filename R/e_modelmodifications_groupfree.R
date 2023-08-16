# Function to free parameters across groups:
groupfree <- function(
  x, # Model
  matrix, # Must be given
  row, # Can be missing
  col, # Can be missing
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
  
  # Obtain rows:
  if (missing(row)){
    row <- unique(x@parameters$row[x@parameters$matrix == matrix])
  }
  if (missing(col)){
    col <- unique(x@parameters$col[x@parameters$matrix == matrix])
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
  
  # If the matrix is symmetric, make sure col contains the lower number:
  if (x@matrices$symmetric[x@matrices$name == matrix]){
    row0 <- row
    col0 <- col
    
    col <- pmin(row0,col0)
    row <- pmax(row0,col0)
    
    # row <- c(row0,col0)
    # col <- c(col0,row0)
  }
  
  # Check if consistent across groups:
  cons <- x@parameters %>% dplyr::group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% dplyr::summarize(consistent = all(.data[['fixed']])|all(!.data[['fixed']]))
  cons <- cons$consistent[cons$matrix %in% matrix & cons$row %in% row & cons$col %in% col]
  
  # which are to be freed:
  whichFree <- which(x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & !x@parameters$fixed)
  
  # Length0?
  if (length(cons)==0 | length(whichFree) == 0){
    if (verbose) message("No parameters need to be freed")
      return(x)
  }
  
  # Any inconsistent?
  if (any(!cons)){
    stop("Model is not consistent across groups. Run unionmodel() or intersectionmodel() first, or manually set parameters consistent across groups.")
  }
  
  # current max par:
  curMax <- max(x@parameters$par)
  
  # Set the model to be not computed:
  x@computed <- FALSE
  
  # Add new parameter values:
  x@parameters$par[whichFree] <- 
    curMax + seq_len(length(whichFree))

  # Relabel:
  x@parameters <- parRelabel(x@parameters)
  

  

  # Clear the parameters (let's keep the old estimates, why not):
  # x@parameters$est[whichFree] <- 0
  # x@parameters$std[whichFree] <- 0
  # x@parameters$se[whichFree] <- 0
  # x@parameters$p[whichFree] <- 0
  # x@parameters$mi[whichFree] <- NA
  # x@parameters$pmi[whichFree] <- NA
  # x@parameters$mi_equal[whichFree] <- NA
  # x@parameters$pmi_equal[whichFree] <- NA
  x@parameters <- clearpars(x@parameters,whichFree)

  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  
  # Output:
  if (verbose){
    message(paste0("Freed ",max(x@parameters$par)-curMax," parameters!"))
  }
  
  
  # Write to log:
  if (log){
    # Add log:
    x <- addLog(x, paste0("Set element(s) of ",matrix," free across groups. Freed ",max(x@parameters$par) - curMax," parameters!")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  

  return(x)
}