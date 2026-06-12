setestimator <- function(x, estimator){
  # ml_lvm models cannot be switched to/from the two-level sufficient
  # statistics ML estimator after creation, as this changes the distribution
  # ("TwoLevelGaussian") and requires (re)computing the two-level sufficient
  # statistics. Rebuild the model instead:
  if (x@model == "ml_lvm" && estimator != x@estimator &&
      (estimator == "ML" || (length(x@distribution) == 1 && x@distribution == "TwoLevelGaussian"))){
    stop("The estimator of an ml_lvm model cannot be switched between 'FIML' and the two-level 'ML' estimator after the model is created. Rebuild the model with ml_lvm(..., estimator = '", estimator, "').")
  }
  # Validate that the sample slots required by the new estimator are present:
  .check_estimator_requirements <- function(model, estimator){
    if (estimator == "FIML"){
      if (length(model@sample@fimldata) == 0){
        stop("estimator = 'FIML' requires raw data (per-observation FIML data). ",
             "The model was not built from raw data, so FIML cannot be used. ",
             "Rebuild the model from a raw data frame (the `data` argument) instead of from covariances.")
      }
    } else if (estimator %in% c("ULS","WLS","DWLS")){
      if (length(model@sample@WLS.W) == 0){
        stop("estimator = '", estimator, "' requires a weights matrix (WLS.W) in the sample. ",
             "This is only available when the model was built with WLS.W/Gamma sample statistics ",
             "(e.g. from raw data with ordered/ordinal variables or by supplying WLS.W). ",
             "Rebuild the model accordingly before switching to '", estimator, "'.")
      }
    }
    invisible(TRUE)
  }

  # Recursively set the estimator on a model and its baseline/saturated:
  .set_one <- function(model, estimator){
    .check_estimator_requirements(model, estimator)
    model@estimator <- estimator
    model@computed <- FALSE
    # Propagate to baseline and saturated (themselves psychonetrics models):
    if (!is.null(model@baseline_saturated$baseline) &&
        is(model@baseline_saturated$baseline, "psychonetrics")){
      model@baseline_saturated$baseline <- .set_one(model@baseline_saturated$baseline, estimator)
    }
    if (!is.null(model@baseline_saturated$saturated) &&
        is(model@baseline_saturated$saturated, "psychonetrics")){
      model@baseline_saturated$saturated <- .set_one(model@baseline_saturated$saturated, estimator)
    }
    model
  }

  x <- .set_one(x, estimator)
  x
}
