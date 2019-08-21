getVCOV <- function(model){
  if (!is.null(model@information)){
    Info <- model@information
  } else {
    Info <- psychonetrics_FisherInformation(model)
  }
  
  1/sum(model@sample@groups$nobs) * as.matrix(solve_symmetric(Info))
}