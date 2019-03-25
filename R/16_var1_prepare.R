# Prepare all matrices for the fit, gradient and hessian of var1 models:
prepare_var1 <- function(x, model){

  # New model:
  newMod <- updateModel(x,model)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Form the model matrices:
  # mats <- formModelMatrices(newMod)
  
  # Add a few that will be useful:
  # for (g in 1:nGroup){
  #   IminO <- Diagonal(nrow(mats[[g]]$omega_zeta)) - mats[[g]]$omega_zeta
  # 
  #   if (any(eigen(IminO)$values < 0)){
  #     warning("I - Omega_zeta is not positive definite, gradient may not be correct.")
  #     mats[[g]]$OmegaStar <- solve_symmetric(IminO)
  #   }  else {
  #     mats[[g]]$OmegaStar <- solve_symmetric(IminO)
  #   }
  #   mats[[g]]$DeltaOmegaStar <- mats[[g]]$delta_zeta %*% mats[[g]]$OmegaStar 
  #   mats[[g]]$BetaStar <- solve_symmetric(Diagonal(nrow(mats[[g]]$beta)^2) - (mats[[g]]$beta %x% mats[[g]]$beta))
  #   mats[[g]]$L_betaStar <- model@extramatrices$L %*%  mats[[g]]$BetaStar 
  #   mats[[g]]$IkronBeta <- model@extramatrices$In %x% mats[[g]]$beta
  #   mats[[g]]$E <- Emat(nrow(mats[[g]]$beta),mats[[g]]$beta)
  #   mats[[g]]$SigmaZeta <- as.matrix(mats[[g]]$delta_zeta %*% mats[[g]]$OmegaStar %*% mats[[g]]$delta_zeta)
  #   mats[[g]]$sigmaZetaVec <- as.vector(mats[[g]]$SigmaZeta)
  #   if (model@rawts){
  #     mats[[g]]$mis <-  model@sample@missingness[[g]]  
  #   }
  #   
  # }
  
  # Raw ts?
  # if (!model@rawts){
    imp <- implied_var1(newMod,all=FALSE)
  # } else {
  #   stop("rawts is not supported")
  #   # Compute implied matrices:
  #   imp <- implied_var1_rawts(newMod)
  # }
  
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
    # if (model@rawts){
      # mats[[g]] <- mats[[g]][names(mats[[g]]) != "mu"]
    # }
    groupModels[[g]] <- c( imp[[g]], mMat, model@extramatrices, model@types) # FIXME: This will lead to extra matrices to be stored?
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