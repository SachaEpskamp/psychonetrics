# Function to prune all non-significant results:
prune <- function(
  x, # Model
  alpha = 0.01, # Significance
  bonferroni = FALSE,
  matrices, # Automatically chosen
  runmodel = TRUE,
  recursive = TRUE,
  verbose = TRUE,
  log = TRUE,
  ...){
  # If not computed, nothing to do:
  if (!x@computed){
    stop("Model must have been computed first.")
  }
  
  
  if (!runmodel & recursive){
    stop("recursive = TRUE requires runmodel = TRUE")
  }
  # Select the matrices:
  if (missing(matrices)){
    if (x@model == "ggm"){
      matrices <- "omega"
    } else if (x@model == "precision"){
      matrices <- "kappa"
    }
  }
  
  # Which parameters to test:
  whichTest <- which(x@parameters$matrix %in% matrices & !x@parameters$fixed & x@parameters$var1_id!=x@parameters$var2_id)
  
  # Number of tests:
  nTest <- length(unique(x@parameters$par[whichTest]))

  # If no tests, break:
  if (nTest == 0){
    return(x)
  }
  
  # Adjust alpha:
  if (bonferroni){
    alpha_adjust <- alpha / nTest    
  } else {
    alpha_adjust <- alpha
  }
  
  # Test for significance:
  nonsig <- x@parameters$p > alpha_adjust
  
  # If any non-sig, adjust:
  if (all(is.na(nonsig)) || !any(nonsig[!is.na(nonsig)])){
    if (log){
      # Add log:
      x <- addLog(x, paste("Pruned all parameters in matrices ",paste0(matrices,collapse=", ")," at alpha = ",alpha," (none were removed)")) 
    }    
    return(x)
  } 
  
  if (verbose){
    message(paste0("Clearing ",length(unique(x@parameters$par[nonsig]))," parameters!"))
  }
  
  # Set computed:
  x@computed <- FALSE
  
  # Clear the parameters:
  x@parameters$est[nonsig] <- 0
  x@parameters$std[nonsig] <- 0
  x@parameters$se[nonsig] <- 0
  x@parameters$fixed[nonsig] <- TRUE
  x@parameters$par[nonsig] <- 0
  x@parameters$mi[nonsig] <- NA
  x@parameters$pmi[nonsig] <- NA
  x@parameters$mi_equal[nonsig] <- NA
  x@parameters$pmi_equal[nonsig] <- NA

  x@parameters   <- parRelabel(x@parameters)
  
  if (log){
    # Add log:
    x <- addLog(x, paste("Pruned all parameters in matrices ",paste0(matrices,collapse=", ")," at alpha = ",alpha_adjust)) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }
  
  # Recurse if needed:
  if (recursive){
    x <- x %>% prune(
      alpha = alpha,
      bonferroni = bonferroni,
      matrices = matrices, # Automatically chosen
      runmodel = TRUE,
      recursive = TRUE,
      verbose = verbose,
      ...
    )
  }
  

  x
}