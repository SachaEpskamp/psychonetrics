# General log-likelihood!
psychonetrics_logLikelihood <- function(model){
  # The two-level ML estimator only has an R implementation:
  model <- force_R_path_if_needed(model)

  # What distribution to use:
  distribution <- model@distribution
  
  # Prepare model:
  prep <- prepareModel(parVector(model), model)
  
  
  # Obtain log likelihood:
  loglikFun <- switch(
    distribution,
    "Gaussian" = logLikelihood_gaussian,
    "Spin" = logLikelihood_Ising,
    "TwoLevelGaussian" = logLikelihood_gaussian2L
  )
  
  # Run and return:
  loglikFun(prep)
}