# General log-likelihood!
psychonetrics_logLikelihood <- function(model){
  # What distribution to use:
  distribution <- model@distribution
  
  # Obtain log likelihood:
  loglikFun <- switch(
    distribution,
    "Gaussian" = logLikelihood_gaussian
  )
  
  # Run and return:
  loglikFun(model)
}