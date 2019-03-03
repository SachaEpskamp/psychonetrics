# Prepare all matrices for the fit, gradient and hessian of gvar models:
prepare_gvar <- function(x, model){

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
    mats[[g]]$OmegaStar <- solve(Diagonal(nrow(mats[[g]]$omega_zeta)) - mats[[g]]$omega_zeta)
    mats[[g]]$DeltaOmegaStar <- mats[[g]]$delta_zeta %*% mats[[g]]$OmegaStar 
    mats[[g]]$BetaStar <- solve(Diagonal(nrow(mats[[g]]$beta)^2) - (mats[[g]]$beta %(x)% mats[[g]]$beta))
    mats[[g]]$L_betaStar <- model@extramatrices$L %*%  mats[[g]]$BetaStar 
    mats[[g]]$IkronBeta <- model@extramatrices$In %(x)% mats[[g]]$beta
    mats[[g]]$E <- Emat(nrow(mats[[g]]$beta),mats[[g]]$beta)
    mats[[g]]$SigmaZeta <- as.matrix(mats[[g]]$delta_zeta %*% mats[[g]]$OmegaStar %*% mats[[g]]$delta_zeta)
    mats[[g]]$sigmaZetaVec <- as.vector(mats[[g]]$SigmaZeta)
    if (model@rawts){
      mats[[g]]$mis <-  model@sample@missingness[[g]]  
    }
    
  }
  
  # Raw ts?
  if (!model@rawts){
    imp <- implied_gvar_norawts(mats)
  } else {
    # Compute implied matrices:
    imp <- implied_gvar_rawts(mats)
  }
  
  # Sample stats:
  S <- model@sample@covs
  means <- model@sample@means
  nVar <- nrow(model@sample@variables)
  
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
    if (model@rawts){
      mats[[g]] <- mats[[g]][names(mats[[g]]) != "mu"]
    }
    groupModels[[g]] <- c( mats[[g]],imp[[g]], mMat, model@extramatrices) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    if (model@rawts){
      groupModels[[g]]$Drawts <- model@Drawts[[g]]
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