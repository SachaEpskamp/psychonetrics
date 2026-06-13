usecpp <- function(x, use = TRUE){
  # The two-level ML estimator (estimator = "ML") with within-cluster missing
  # data is only implemented on the R path (the per-pattern missing-data
  # likelihood has no C++ twin). Enabling C++ for such a model would call a
  # C++ prepare that lacks the missing-data structures and fail cryptically, so
  # we keep it on the R path and inform the user instead:
  if (isTRUE(use) && twolevel_model_has_missing(x)){
    message("Note: the two-level ML estimator with missing data is evaluated ",
            "in R only; keeping cpp = FALSE for this model.")
    use <- FALSE
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