# Robust ML estimators (Phase 1: complete data only). Each maps internally to
# estimator = "ML" (point estimates are plain ML, identical to estimator="ML")
# plus a robust configuration controlling the standard errors and the scaled
# test statistic. See get_robust_config() in 00_codeOrganization.R.
#   MLM   : Browne (1984) robust.sem SEs + Satorra-Bentler (1994) scaled chi-square
#   MLMV  : robust.sem SEs + scaled-and-shifted chi-square (Satorra 2000; Asparouhov & Muthen 2010)
#   MLMVS : robust.sem SEs + mean-and-variance-adjusted (Satterthwaite) chi-square
#   MLR   : Huber-White (sandwich) SEs + Yuan-Bentler-Mplus scaled chi-square
.robustEstimators <- c("MLM","MLMV","MLMVS","MLR")

.robust_config_for <- function(estimator){
  switch(estimator,
    "MLM"   = list(se = "robust.sem",         test = "satorra.bentler",    label = "MLM"),
    "MLMV"  = list(se = "robust.sem",         test = "scaled.shifted",     label = "MLMV"),
    "MLMVS" = list(se = "robust.sem",         test = "mean.var.adjusted",  label = "MLMVS"),
    "MLR"   = list(se = "robust.huber.white", test = "yuan.bentler.mplus", label = "MLR")
  )
}

# Resolve a (possibly robust) estimator name supplied to a model constructor
# into the internal estimator and the robust configuration. Used by lvm() /
# varcov() (and the families that wrap them) so that estimator = "MLM" etc. is
# accepted directly at construction, equivalent to building with
# estimator = "ML" and then calling setestimator(model, "MLM").
# Returns list(estimator = <internal>, robust = <config or empty list>).
resolve_robust_estimator <- function(estimator){
  if (length(estimator) == 1 && estimator %in% .robustEstimators){
    list(estimator = "ML", robust = .robust_config_for(estimator))
  } else {
    list(estimator = estimator, robust = list())
  }
}

setestimator <- function(x, estimator){
  # ml_lvm models cannot be switched to/from the two-level sufficient
  # statistics ML estimator after creation, as this changes the distribution
  # ("TwoLevelGaussian") and requires (re)computing the two-level sufficient
  # statistics. Rebuild the model instead:
  if (x@model == "ml_lvm" && estimator != x@estimator &&
      (estimator == "ML" || (length(x@distribution) == 1 && x@distribution == "TwoLevelGaussian"))){
    stop("The estimator of an ml_lvm model cannot be switched between 'FIML' and the two-level 'ML' estimator after the model is created. Rebuild the model with ml_lvm(..., estimator = '", estimator, "').")
  }

  # Robust ML estimators map internally to estimator = "ML" plus a robust
  # configuration. Resolve the requested name into the internal estimator and
  # the robust config now:
  if (estimator %in% .robustEstimators){
    robust_cfg <- .robust_config_for(estimator)
    internal_estimator <- "ML"
  } else {
    robust_cfg <- list()
    internal_estimator <- estimator
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

  # Robust ML (Phase 1) requires raw data: the sandwich corrections need the
  # asymptotic covariance Gamma (MLM family) or casewise scores (MLR), neither
  # of which is available from summary statistics alone. FIML (missing data)
  # support is deferred to Phase 2.
  if (length(robust_cfg) > 0){
    has_raw <- nrow(x@sample@rawdata) > 0 || length(x@sample@fimldata) > 0
    if (!has_raw){
      stop("estimator = '", robust_cfg$label, "' (robust ML) requires raw data, which ",
           "is not stored in this model. Rebuild the model from a raw data frame ",
           "(the `data` argument) before switching to '", robust_cfg$label, "'.")
    }
    if (length(x@sample@fimldata) > 0 && nrow(x@sample@rawdata) == 0){
      stop("estimator = '", robust_cfg$label, "' currently supports complete data only ",
           "(Phase 1). The model contains missing data (FIML). Robust ML with FIML is ",
           "planned for a future release.")
    }
  }

  # Recursively set the estimator on a model and its baseline/saturated:
  .set_one <- function(model, estimator, robust_cfg){
    .check_estimator_requirements(model, estimator)
    model@estimator <- estimator
    if (.hasSlot(model, "robust")){
      model@robust <- robust_cfg
    }
    model@computed <- FALSE
    # Propagate to baseline and saturated (themselves psychonetrics models):
    if (!is.null(model@baseline_saturated$baseline) &&
        is(model@baseline_saturated$baseline, "psychonetrics")){
      model@baseline_saturated$baseline <- .set_one(model@baseline_saturated$baseline, estimator, robust_cfg)
    }
    if (!is.null(model@baseline_saturated$saturated) &&
        is(model@baseline_saturated$saturated, "psychonetrics")){
      model@baseline_saturated$saturated <- .set_one(model@baseline_saturated$saturated, estimator, robust_cfg)
    }
    model
  }

  x <- .set_one(x, internal_estimator, robust_cfg)
  x
}
