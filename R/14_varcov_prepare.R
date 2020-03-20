# Prepare all matrices for the fit, gradient and hessian of varcov models:
prepare_varcov <- function(x, model){
  # New model:
  if (model@cpp){
    newMod <- updateModel_cpp(x,model,FALSE)
  } else {
    newMod <- updateModel(x,model,updateMatrices = FALSE)  
  }
  
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Form the model matrices:
  # mats <- formModelMatrices(newMod)
  
  # Compute implied matrices:
  if (model@cpp){
    imp <- implied_varcov_cpp(newMod, all = FALSE)  
  } else {
    imp <- implied_varcov(newMod, all = FALSE)
  }
  

  # Sample stats:
  S <- model@sample@covs
  means <- model@sample@means
  nVar <- nrow(model@sample@variables)
  thresholds <- model@sample@thresholds
  
  # Extra mats:
  mMat <- list(
    M = Mmatrix(model@parameters)
  )

 # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    groupModels[[g]] <- c(imp[[g]], mMat, model@extramatrices, model@types) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    groupModels[[g]]$corinput <- model@sample@corinput
    groupModels[[g]]$meanstructure <- model@meanstructure
    groupModels[[g]]$cpp <- model@cpp
    
    if (is.null( groupModels[[g]]$tau)){
      groupModels[[g]]$tau <- matrix(NA,1,nVar)
    }
  
    if (length(thresholds) > 0){
      groupModels[[g]]$thresholds <- thresholds[[g]]      
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