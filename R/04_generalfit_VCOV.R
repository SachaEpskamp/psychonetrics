getVCOV <- function(model,approximate_SEs = FALSE){
  if (!is.null(model@information)){
    Info <- model@information
  } else {
    if (model@cpp){
      Info <- psychonetrics_FisherInformation_cpp(model)
    } else {
      Info <- psychonetrics_FisherInformation(model)      
    }

  }
  
  # Attempt to invert:
  1/sum(model@sample@groups$nobs) * as.matrix(solve_symmetric(Info,approx=approximate_SEs))
}