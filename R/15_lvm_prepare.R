# Prepare all matrices for the fit, gradient and hessian of lvm models:
prepare_lvm <- function(x, model){

  # New model:
  newMod <- updateModel(x,model,updateMatrices = FALSE)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Form the model matrices:
  mats <- formModelMatrices(newMod)
  
  # Compute implied matrices:
  imp <- implied_lvm(newMod, all = FALSE)
 
  # Sample stats:
  S <- model@sample@covs
  means <- model@sample@means
  nVar <- nrow(model@sample@variables)
  thresholds <- model@sample@thresholds
  
  # Extra mats:
  mMat <- list(
    M = Mmatrix(model@parameters)
  )
  # extraMatrices <- list(
  #   M = Mmatrix(model@parameters), # Model matrix
  #   D = duplicationMatrix(nVar), # non-strict duplciation matrix
  #   L = eliminationMatrix(nVar), # Elinimation matrix
  #   Dstar = duplicationMatrix(nVar,diag = FALSE), # Strict duplicaton matrix
  #   A = diagonalizationMatrix(nVar), # Diagonalization matrix
  #   An2 = diagonalizationMatrix(nVar^2), # Larger diagonalization matrix
  #   In = Diagonal(nVar), # Identity of dim n
  #   In2 = Diagonal(nVar^2), # Identity of dim n^2
  #   In3 = Diagonal(nVar^3), # Identity of dim n^3
  #   E = basisMatrix(nVar) # Basis matrix
  # )

  # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    groupModels[[g]] <- c(imp[[g]], mMat, model@extramatrices, model@types) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    # groupModels[[g]]$thresholds <- thresholds[[g]]
  }
  
  
  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal=nTotal,
    nGroup=nGroup,
    groupModels=groupModels
  ))
}