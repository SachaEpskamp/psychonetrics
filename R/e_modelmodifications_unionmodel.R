# Function to prune all non-significant results:
unionmodel <- function(
  x, # Model
  runmodel = FALSE,
  verbose = TRUE,
  log = TRUE,
  identify = TRUE,
  ...){
  # If only one group, nothing to do!
  if (nrow(x@sample@groups) == 1){
    if (verbose) message("Only one group in model, unionmodel does nothing...")
    return(x)
  }
  
  # Copy parameter table:
  pars <- x@parameters
  
  # Add dummy id:
  pars$id <- seq_len(nrow(pars))
  
  # Find set of parameters at least in one group:
  pars <- pars %>%  left_join(pars %>% dplyr::group_by_("matrix","row","col") %>% dplyr::summarise_(anyFree = ~any(!fixed)),
                    by = c("matrix","row","col"))
  
  # So which to free?
  whichFree <- which(pars$fixed & pars$anyFree)
  
  # If nothing to do, return
  if (length(whichFree)==0){
    return(x)
  }
  

  # Set computed:
  x@computed <- FALSE
  
  # Free the parameters:
  x@parameters$fixed[whichFree] <- FALSE
  x@parameters$par[whichFree] <- max(x@parameters$par) + seq_len(length(whichFree))
  # x@parameters$mi[whichFree] <- NA
  # x@parameters$pmi[whichFree] <- NA
  # x@parameters$mi_equal[whichFree] <- NA
  # x@parameters$pmi_equal[whichFree] <- NA
  
  x@parameters <- clearpars(x@parameters,whichFree)

  # Just in case... # FIXME: Not needed right?
  x@parameters   <- parRelabel(x@parameters)
  
  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  if (verbose){
    message(paste0("Freeing ",length(whichFree)," parameters!"))
  }
  
  
  if (log){
    # Add log:
    x <- addLog(x, paste("Unified models across groups (union model). Freed in total ",length(whichFree)," parameters")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }

  x
}