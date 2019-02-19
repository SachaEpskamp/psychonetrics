# Stepwise up search:
stepup <- function(
  x, # psychonetrics model
  alpha = 0.01, # Alpha to use for modification indices
  matrices, # Matrices to search
  ... # Fit arguments
){

    if (missing(matrices)){
      if (x@model == "ggm"){
        matrices <- "omega"
      } else if (x@model == "precision"){
        matrices <- "kappa"
      } else stop("No default argument for 'matrices' for current model.")
    }
  
  # Start loop:
  repeat{
  # Stepwise up?
    if (!"network" %in% x@equal){ # FIXME:  I am lazy, must be more general
      if (any(x@parameters$mi[x@parameters$matrix %in% matrices] > qchisq(alpha,1,lower.tail=FALSE))){

        # FIXME: Make nice free parameter function
        best <- which(x@parameters$mi == max(x@parameters$mi[x@parameters$matrix %in% matrices]))
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
      if (any(x@parameters$mi_equal[x@parameters$matrix %in% matrices] > qchisq(alpha,1,lower.tail=FALSE))){
        # FIXME: Make nice free parameter function
        best <- which(x@parameters$mi_equal == max(x@parameters$mi_equal[x@parameters$matrix %in% matrices]))
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
  
    
  }
  
  # Add log:
  x <- addLog(x, "Performed step-up model search")
  
  return(x)
    
}