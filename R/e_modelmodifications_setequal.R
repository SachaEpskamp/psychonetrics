# Function to fix a certain parameter
parequal <- function(
  x, # Model
  ...,
 inds = integer(0), # Indices to set equal
  verbose,
  log = TRUE,
  runmodel = FALSE){
  
  if (missing(verbose)){
    verbose <- x@verbose
  }

  # Capture parameter indices from ...: only unnamed elements are treated as
  # indices to constrain equal. Named elements (e.g. runmodel arguments such
  # as addMIs) are ignored here so they do not leak into the constraint set.
  dots <- list(...)
  if (length(dots) > 0){
    nms <- names(dots)
    if (is.null(nms)){
      dotInds <- unlist(dots)
    } else {
      dotInds <- unlist(dots[nms == ""])
    }
  } else {
    dotInds <- NULL
  }
  inds <- c(dotInds, inds)
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
    # Do NOT forward ... to runmodel: ... is reserved for parameter indices,
    # and numeric indices would otherwise leak into runmodel positionally.
    x <- x %>% runmodel(verbose=verbose)
  }
  

  return(x)
}