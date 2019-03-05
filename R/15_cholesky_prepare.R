# Prepare all matrices for the fit, gradient and hessian of cholesky models:
prepare_cholesky <- function(x, model){
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
  
  # Compute implied matrices:
  imp <- implied_cholesky(mats)

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
  }
  
  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal=nTotal,
    nGroup=nGroup,
    groupModels=groupModels
  ))
}