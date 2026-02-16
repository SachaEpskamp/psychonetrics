# Prepare all matrices for the fit, gradient and hessian of meta_var1 models:
prepare_meta_var1 <- function(x, model){
  # New model:
  newMod <- updateModel(x, model, updateMatrices = FALSE)

  # Number of groups:
  nGroup <- nrow(model@sample@groups)

  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)

  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Form the model matrices:
  mats <- formModelMatrices(newMod)

  # Compute implied matrices:
  imp <- implied_meta_var1(newMod, all = FALSE)

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
    groupModels[[g]] <- c(imp[[g]], mMat, model@extramatrices, model@types)
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    groupModels[[g]]$corinput <- FALSE # The meta-level data is not corinput
    groupModels[[g]]$cpp <- model@cpp
    groupModels[[g]]$meanstructure <- model@meanstructure
    # Override D with D_c: the ML gradient uses D for the meta-level sigma (nCov x nCov),
    # not the VAR1-level D (nNode x nNode):
    groupModels[[g]]$D <- model@extramatrices$D_c
  }

  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal = nTotal,
    nGroup = nGroup,
    groupModels = groupModels
  ))
}
