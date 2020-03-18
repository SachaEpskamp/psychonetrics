# General fit function!
psychonetrics_fitfunction <- function(x, model){
  
  # Prepare model:
  prep <- prepareModel(x, model)
  
  # What estimator to use:
  estimator <- model@estimator
  
  # What distribution to use:
  distribution <- model@distribution

  # Obtain estimator:
  if (model@cpp){
    estFun <- switch(
      estimator,
      "ML" = switch(distribution,"Gaussian" = maxLikEstimator_Gauss, "Ising" = maxLikEstimator_Ising),
      "ULS" =  switch(distribution,"Gaussian" = ULS_Gauss),
      "DWLS" = switch(distribution,"Gaussian" = ULS_Gauss),
      "WLS" = switch(distribution,"Gaussian" = ULS_Gauss),
      "FIML" = switch(distribution,"Gaussian" = fimlEstimator_Gauss)
    )
  } else {
    estFun <- switch(
      estimator,
      "ML" = switch(distribution,"Gaussian" = maxLikEstimator_Gauss, "Ising" = maxLikEstimator_Ising),
      "ULS" =  switch(distribution,"Gaussian" = ULS_Gauss),
      "DWLS" = switch(distribution,"Gaussian" = ULS_Gauss),
      "WLS" = switch(distribution,"Gaussian" = ULS_Gauss),
      "FIML" = switch(distribution,"Gaussian" = fimlEstimator_Gauss)
    )
  }
  # Run and return:
  estFun(prep)
}