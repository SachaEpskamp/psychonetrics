# Prepare all matrices for the fit, gradient and hessian of ggm models:
prepare_ggm <- function(x, model){
  # New model:
  newMod <- updateModel(x,model)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Form the model matrices:
  mats <- formModelMatrices(newMod)
  
  # Add a few that will be useful:
  for (g in 1:nGroup){
    mats[[g]]$IminOinv <- solve(Diagonal(nrow(mats[[g]]$omega)) - mats[[g]]$omega)
  }

  # Compute implied matrices:
  imp <- implied_ggm(mats)

  
  # Sample stats:
  S <- model@sample@covs
  means <- model@sample@means
  nVar <- nrow(model@sample@variables)
  
  # Extra mats:
  mMat <- list(
    M = Mmatrix(model@parameters)
  )
   # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    groupModels[[g]] <- c(mats[[g]],imp[[g]], mMat, model@extramatrices) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    if (model@rawts){
      groupModels[[g]]$Drawts <- model@Drawts[[g]]
    }
  }
  
  # rawts
  if (model@rawts){
    # Change the implied matrices to full size:
    for (g in 1:nGroup){
      # Missing pattern:
      missings <- model@sample@missingness[[g]]
      
      # Create the massive matrix:
      muFull <- Reduce("rbind",lapply(seq_len(nrow(missings)),function(x)groupModels[[g]]$mu))
      sigFull <- Reduce("bdiag",lapply(seq_len(nrow(missings)),function(x)groupModels[[g]]$sigma))
  
      # Cut out the rows and cols:
      obsvec <- !as.vector(t(missings))
      muFull <- muFull[obsvec,,drop=FALSE]
      sigFull <- sigFull[obsvec,obsvec]
      
      # Overwrite:
      groupModels[[g]]$mu <- muFull
      groupModels[[g]]$sigma <- sigFull
      groupModels[[g]]$kappa <- solve(sigFull)
      # groupModels[[g]]$missings <- missings
      # groupModels[[g]]$rawts <- TRUE
      # groupModels[[g]]$blockdiagonal <- TRUE
    }
    
  }
  
  
  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal=nTotal,
    nGroup=nGroup,
    groupModels=groupModels
  ))
}