# Prepare all matrices for the fit, gradient and hessian of panelvar models:
prepare_panelvar <- function(x, model){

  # New model:
  newMod <- updateModel(x,model)

  # Number of groups:
  nGroup <- nrow(model@sample@groups)

  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)

  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Implied model:
  imp <- implied_panelvar(newMod,all=FALSE)

  # Sample stats:
  S <- model@sample@covs
  means <- model@sample@means

  # Extra mats:
  mMat <- list(
    M = Mmatrix(model@parameters)
  )

  # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    groupModels[[g]] <- c( imp[[g]], mMat, model@extramatrices, model@types)
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
