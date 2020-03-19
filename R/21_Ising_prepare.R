# Prepare all matrices for the fit, gradient and hessian of Ising models:
prepare_Ising <- function(x, model){
  
  # New model:
  newMod <- updateModel(x,model,updateMatrices = FALSE)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs
  
  # Form the model matrices:
  # mats <- formModelMatrices(newMod)
  
  # Compute implied matrices (not needed here, implied is same as mats):
  imp <- implied_Ising(newMod, all = FALSE)
  # imp <- mats
  
  # Sample stats:
  squares <- model@sample@squares
  means <- model@sample@means
  nobs <- model@sample@groups$nobs
  nVar <- nrow(model@sample@variables)
  
  # Extra mats:
  mMat <- list(
    M = Mmatrix(model@parameters)
  )
  
  # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    groupModels[[g]] <- c(imp[[g]], mMat, model@extramatrices) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$squares <- squares[[g]]
    groupModels[[g]]$means <- means[[g]]
    groupModels[[g]]$nobs <- nobs[[g]]
    groupModels[[g]]$cpp <- model@cpp
    
    # Expectations:
    exp <-  isingExpectation(
      groupModels[[g]]$omega,
      groupModels[[g]]$tau,
      groupModels[[g]]$beta,
      groupModels[[g]]$responses
    )    
    groupModels[[g]]$Z <- exp$Z
    groupModels[[g]]$exp_v1 <- exp$exp_v1
    groupModels[[g]]$exp_v2 <- exp$exp_v2
    groupModels[[g]]$exp_H <- exp$exp_H
   
  }
  
  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal=nTotal,
    nGroup=nGroup,
    groupModels=groupModels
  ))
}