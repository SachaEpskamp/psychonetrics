# Function to fix parameters equal across groups:
groupequal <- function(
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
  cons <- x@parameters %>% dplyr::group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% dplyr::summarize(consistent = all(.data[['fixed']])|all(!.data[['fixed']]), anyIdentified = any(.data[['identified']]))
  cons <- cons$consistent[cons$matrix %in% matrix & cons$row %in% row & cons$col %in% col & !cons$anyIdentified]
  
  # which are constrained:
  whichCons <- which(x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & !x@parameters$fixed)
  
  # Length?
  if (length(whichCons) == 0){
    if (verbose) message("No parameters need to be constrained...")
      return(x)
  }
  
  # Any inconsistent?
  if (length(cons) > 0 && any(!cons)){
    stop("Model is not consistent across groups. Run unionmodel() or intersectionmodel() first, or manually set parameters consistent across groups.")
  }
  
  # current max par:
  curMax <- max(x@parameters$par)
  
  # Set the model to be not computed:
  x@computed <- FALSE
  
  # for each group, set parameters equal to first group FIXME: Not pretty but gets the job done:
  id1 <- x@sample@groups$id[1]
  for (g in x@sample@groups$id[-1]){
    group_g <- x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & x@parameters$group_id == g
    group_1 <- x@parameters$matrix == matrix & x@parameters$row %in% row & x@parameters$col %in% col & x@parameters$group_id == id1

    x@parameters$par[group_g] <- x@parameters$par[group_1]
    x@parameters$est[group_g] <- x@parameters$est[group_1]
    x@parameters$fixed[group_g] <- x@parameters$fixed[group_1]
    # x@parameters$identified[group_g] <- x@parameters$identified[group_1]
  }
  

  # Relabel:
  x@parameters   <- parRelabel(x@parameters)
  
  # Clear the parameters:
  # x@parameters$est[whichCons] <- 0
  # x@parameters$std[whichCons] <- 0
  # x@parameters$se[whichCons] <- NA
  # x@parameters$p[whichCons] <- NA
  # x@parameters$mi[whichCons] <- NA
  # x@parameters$pmi[whichCons] <- NA
  # x@parameters$mi_equal[whichCons] <- NA
  # x@parameters$pmi_equal[whichCons] <- NA
  x@parameters <- clearpars(x@parameters,whichCons)
  
  # For every parameter, take mean as starting value:
  x@parameters  <- x@parameters %>% group_by(.data[["par"]]) %>% mutate(est=ifelse(.data[['par']]==0,.data[['est']],mean(.data[['est']],na.rm=TRUE))) %>% ungroup

  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  # Output:
  if (verbose){
    message(paste0("Constrained ",curMax - max(x@parameters$par)," parameters!"))
  }
  
  
  # Write to log:
  if (log){
    # Add log:
    x <- addLog(x, paste0("Set element(s) of ",matrix," equal across groups. Constrained ",curMax - max(x@parameters$par)," parameters!")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  

  return(x)
}