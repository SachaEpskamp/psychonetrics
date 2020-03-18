# maxLikEstimator <- function(x, model){
#   # What distribution?
#   distribution <- model@distribution
#   
#   # Function to use:
#   distFun <- switch(distribution,
#           "Gaussian" = maxLikEstimator_Gauss,
#           "Ising" = maxLikEstimator_Ising)
#   
#   # Run and return:
#   distFun(x, model)
# }