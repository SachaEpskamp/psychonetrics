usecpp <- function(x, use = TRUE){
  # The two-level ML estimator (ml_lvm with estimator = "ML") only has an R
  # implementation; C++ cannot be enabled for these models:
  if (isTRUE(use) && length(x@distribution) == 1 && x@distribution == "TwoLevelGaussian"){
    message("The two-level ML estimator currently only has an R implementation; 'usecpp' is ignored for this model.")
    return(x)
  }
  x@cpp <- use

  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@cpp <- use
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@cpp <- use
  }
  
  x
}