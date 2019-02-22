# Stepwise up search:
stepup <- function(
  x, # psychonetrics model
  alpha = 0.01, # Alpha to use for modification indices
  criterion = "bic", # Stop when criterion is no longer improved. Can also be none to ignore
  matrices, # Matrices to search
  ... # Fit arguments
){

    if (missing(matrices)){
      if (x@model == "ggm"){
        matrices <- "omega"
      } else if (x@model == "precision"){
        matrices <- "kappa"
      }  else if (x@model == "lnm"){
        matrices <- "omega_eta"
      } else stop("No default argument for 'matrices' for current model.")
    }
  
  # Check if MIs are added:
  if (all(is.na(x@parameters$mi))){
    x <- x %>% addMIs(matrices = matrices)
  }
  
  # Start loop:
  repeat{
    oldMod <- x
  # Stepwise up?
    if (!"network" %in% x@equal){ # FIXME:  I am lazy, must be more general
      if (any(x@parameters$mi[x@parameters$matrix %in% matrices & x@parameters$fixed] > qchisq(alpha,1,lower.tail=FALSE))){

        # FIXME: Make nice free parameter function
        best <- which(x@parameters$mi == max(x@parameters$mi[x@parameters$matrix %in% matrices & x@parameters$fixed]))
        x@parameters$par[best] <- max(x@parameters$par) + 1
        x@parameters$fixed[best] <- FALSE
        
        # Update the model:
        x@fitfunctions$extramatrices$M <- Mmatrix(x@parameters) # FIXME: Make nice function for this
        
        # Run:
        x <- x %>% runmodel(...,log=FALSE)
      } else {
        break
      }
    } else {
      if (any(x@parameters$mi_equal[x@parameters$matrix %in% matrices & x@parameters$fixed] > qchisq(alpha,1,lower.tail=FALSE))){
        # FIXME: Make nice free parameter function
        best <- which(x@parameters$mi_equal == max(x@parameters$mi_equal[x@parameters$matrix %in% matrices & x@parameters$fixed]))
        x@parameters$par[best] <- max(x@parameters$par) + 1
        x@parameters$fixed[best] <- FALSE
        
        # Update the model:
        x@fitfunctions$extramatrices$M <- Mmatrix(x@parameters)
        
        # Run:
        x <- x %>% runmodel(...,log=FALSE)
      } else {
        break
      }
    }
    
    # Check criterion:
    if (criterion != "none"){
      oldCrit <- oldMod@fitmeasures[[criterion]]
      newCrit <- x@fitmeasures[[criterion]]
      if (oldCrit < newCrit){
        x <- oldMod
        break
      }
    }
  }
  
  # Add log:
  x <- addLog(x, "Performed step-up model search")
  
  return(x)
    
}