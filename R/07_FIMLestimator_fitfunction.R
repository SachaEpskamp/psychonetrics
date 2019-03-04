fimlEstimator <- function(x, model){
  # What distribution?
  distribution <- model@distribution
  
  # Function to use:
  distFun <- switch(distribution,
          "Gaussian" = fimlEstimator_Gauss)
  
  # Run and return:
  distFun(x, model)
}