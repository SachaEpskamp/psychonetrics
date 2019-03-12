# Prepare all matrices for the fit, gradient and hessian of varcov models:
prepare_varcov <- function(x, model){
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
  if (model@types$y == "ggm"){
    for (g in 1:nGroup){
      mats[[g]]$IminOinv <- corpcor::pseudoinverse(Diagonal(nrow(mats[[g]]$omega)) - mats[[g]]$omega)
    }
  }

  # Compute implied matrices:
  # imp <- implied_varcov(mats)

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
    groupModels[[g]] <- c(mats[[g]], mMat, model@extramatrices, model@types) # FIXME: This will lead to extra matrices to be stored?
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