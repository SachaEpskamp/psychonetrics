# Prepare all matrices for the fit, gradient and hessian of ml_lvm models:
prepare_ml_lvm <- function(x, model){

  # New model:
  newMod <- updateModel(x,model)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Total sample:
  nTotal <- sum(model@sample@groups$nobs)
  
  # Sample per group:
  nPerGroup <- model@sample@groups$nobs

  # Two-level sufficient-statistics ML estimator?
  twolevel_ML <- model@estimator == "ML"

  # Implied model:
  imp <- implied_ml_lvm(newMod, all = FALSE, twolevel_only = twolevel_ML)

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

  # Two-level sufficient statistics (per group), guarded for objects saved
  # before the 'twolevel' slot existed:
  if (twolevel_ML){
    twolevelStats <- get_twolevel_stats(model@sample)
    if (length(twolevelStats) != nGroup){
      stop("estimator = 'ML' for ml_lvm requires two-level sufficient statistics, which are missing from the model. Rebuild the model with ml_lvm(..., estimator = 'ML').")
    }
  }

  # Fill per group:
  groupModels <- list()
  for (g in 1:nGroup){
    # if (model@rawts){
      # mats[[g]] <- mats[[g]][names(mats[[g]]) != "mu"]
    # }
    groupModels[[g]] <- c( imp[[g]], mMat, model@extramatrices, model@types) # FIXME: This will lead to extra matrices to be stored?
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    if (twolevel_ML){
      groupModels[[g]]$twolevel <- twolevelStats[[g]]
    }
  }


  # Return
  return(list(
    nPerGroup = nPerGroup,
    nTotal=nTotal,
    nGroup=nGroup,
    groupModels=groupModels,
    # Flag so that the model Jacobian dispatch (d_phi_theta_ml_lvm) picks the
    # two-level variant with rows [mu; vech Sigma_W; vech Sigma_B]:
    twolevel_ML = twolevel_ML
  ))
}