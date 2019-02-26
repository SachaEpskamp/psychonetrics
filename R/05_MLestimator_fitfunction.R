maxLikEstimator <- function(x, model){
  # What distribution?
  distribution <- model@distribution
  
  # Function to use:
  distFun <- switch(distribution,
          "Gaussian" = maxLikEstimator_Gauss)
  
  # Run and return:
  distFun(x, model)
}