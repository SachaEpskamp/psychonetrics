# Prepare all matrices for the fit, gradient and Hessian of ml_var1 models.
# The framework is always two-level ML (distribution "TwoLevelGaussian"); the
# estimator layer consumes mu, sigma_within, sigma_between, twolevel and D_y.
prepare_ml_var1 <- function(x, model){

  # New model:
  newMod <- updateModel(x, model)

  # Number of groups:
  nGroup <- nrow(model@sample@groups)

  # Total sample (= total number of clusters across groups):
  nTotal <- sum(model@sample@groups$nobs)

  # Sample per group (= number of subjects J_g):
  nPerGroup <- model@sample@groups$nobs

  # Implied model:
  imp <- implied_ml_var1(newMod, all = FALSE)

  # Sample stats (unused by the two-level estimator, but kept parallel to the
  # other frameworks / generic code):
  S <- model@sample@covs
  means <- model@sample@means

  # Model matrix (maps free parameters to the model-matrix vec's):
  mMat <- list(M = Mmatrix(model@parameters))

  # Two-level sufficient statistics per group:
  twolevelStats <- get_twolevel_stats(model@sample)
  if (length(twolevelStats) != nGroup){
    stop("ml_var1 requires two-level sufficient statistics, which are missing from the model. Rebuild the model with ml_var1(...).")
  }

  # Fill per group:
  groupModels <- list()
  for (g in seq_len(nGroup)){
    groupModels[[g]] <- c(imp[[g]], mMat, model@extramatrices, model@types)
    groupModels[[g]]$S <- S[[g]]
    groupModels[[g]]$means <- means[[g]]
    groupModels[[g]]$twolevel <- twolevelStats[[g]]
  }

  # Return:
  return(list(
    nPerGroup = nPerGroup,
    nTotal = nTotal,
    nGroup = nGroup,
    groupModels = groupModels
  ))
}
