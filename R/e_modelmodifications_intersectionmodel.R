# Function to prune all non-significant results:
intersectionmodel <- function(
  x, # Model
  runmodel = FALSE,
  verbose = TRUE,
  log = TRUE,
  identify = TRUE,
  ...){
  # If only one group, nothing to do!
  if (nrow(x@sample@groups) == 1){
    if (verbose) message("Only one group in model, intersectionmodel does nothing...")
    return(x)
  }
  
  # Copy parameter table:
  pars <- x@parameters
  
  # Add dummy id:
  pars$id <- seq_len(nrow(pars))
  
  # Find set of parameters at least in one group:
  pars <- pars %>%  left_join(pars %>% dplyr::group_by_("matrix","row","col") %>% dplyr::summarise_(anyFixed = ~any(fixed)),
                    by = c("matrix","row","col"))
  
  # So which to fix?
  whichFix <- which(!pars$fixed & pars$anyFixed)
  
  # If nothing to do, return
  if (length(whichFix)==0){
    return(x)
  }
  

  # Set computed:
  x@computed <- FALSE
  
  # Fix the parameters:
  x@parameters$fixed[whichFix] <- TRUE
  x@parameters$par[whichFix] <- 0
  # x@parameters$mi[whichFix] <- NA
  # x@parameters$pmi[whichFix] <- NA
  # x@parameters$mi_equal[whichFix] <- NA
  # x@parameters$pmi_equal[whichFix] <- NA
  # x@parameters$est[whichFix] <- 0
  # x@parameters$std[whichFix] <- 0
  # x@parameters$se[whichFix] <- 0
  x@parameters <- clearpars(x@parameters,whichFix)
  
  # Relabel
  x@parameters <- parRelabel(x@parameters)
  
  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  if (verbose){
    message(paste0("Fixing ",length(whichFix)," parameters!"))
  }
  
  
  
  if (log){
    # Add log:
    x <- addLog(x, paste("Unified models across groups (intersection model). Fixed in total ",length(whichFix)," parameters")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }

  x
}