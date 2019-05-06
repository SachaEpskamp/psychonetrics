# Function to fix a certain parameter
parequal <- function(
  x, # Model
  ...,
 inds = integer(0), # Indices to set equal
  verbose = TRUE,
  log = TRUE,
  runmodel = FALSE){

  inds <- c(unlist(list(...),inds))
  # current max par:
  curMax <- max(x@parameters$par)
  
  # Set the model to be not computed:
  x@computed <- FALSE
  
  # Fix the parameters:
  x@parameters$est[inds] <- x@parameters$est[inds[1]]
  # x@parameters$std[whichCons] <- NA
  x@parameters$par[inds] <- x@parameters$par[inds[1]]
  # x@parameters$se[whichCons] <- NA
  # x@parameters$p[whichCons] <- NA
  # x@parameters$mi[whichCons] <- NA
  # x@parameters$pmi[whichCons] <- NA
  # x@parameters$mi_equal[whichCons] <- NA
  # x@parameters$pmi_equal[whichCons] <- NA
  # x@parameters$fixed[whichCons] <- TRUE
  x@parameters <- clearpars(x@parameters,inds)

  # Relabel:
  x@parameters   <- parRelabel(x@parameters)
  
  # Output:
  if (verbose){
    message(paste0("Constrained ",curMax - max(x@parameters$par)," parameters!"))
  }
  

  # Write to log:
  if (log){
    # Add log:
    x <- addLog(x, paste0("Constrained ",curMax - max(x@parameters$par)," parameters!")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  

  return(x)
}