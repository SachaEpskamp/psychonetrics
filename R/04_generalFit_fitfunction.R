# General fit function!
psychonetrics_fitfunction <- function(x, model){
  # What estimator to use:
  estimator <- model@estimator
  
  # Obtain estimator:
  estFun <- switch(
    estimator,
    "ML" = maxLikEstimator,
    "ULS" = ULSestimator
  )
  
  # Run and return:
  estFun(x, model)
}