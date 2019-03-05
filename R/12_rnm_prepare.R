# Prepare all matrices for the fit, gradient and hessian of rnm models:
prepare_rnm <- function(x, model){

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
    mats[[g]]$OmegaStar <- corpcor::pseudoinverse(Diagonal(nrow(mats[[g]]$omega_epsilon)) - mats[[g]]$omega_epsilon)
    mats[[g]]$BetaStar <- corpcor::pseudoinverse(Diagonal(nrow(mats[[g]]$beta)) - mats[[g]]$beta)
    mats[[g]]$Lambda_BetaStar <- mats[[g]]$lambda %*%  mats[[g]]$BetaStar 
    mats[[g]]$Betasta_sigmaZeta <- mats[[g]]$BetaStar %*% mats[[g]]$sigma_zeta
    mats[[g]]$sigma_epsilon <- mats[[g]]$delta_epsilon %*% mats[[g]]$OmegaStar %*% mats[[g]]$delta_epsilon 
    mats[[g]]$tBetakronBeta <- t(mats[[g]]$BetaStar) %(x)% mats[[g]]$BetaStar
  }

  # Compute implied matrices:
  imp <- implied_rnm(mats)

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