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

  # ml_var1 supports only the two-level 'ML' estimator (v1):
  if (x@model == "ml_var1" && estimator != "ML"){
    stop("ml_var1 models only support estimator = 'ML' (the two-level pseudo-ML estimator). For full-information ML use panelvar()/dlvm1().")
  }

  # ml_varcov supports 'ML' (two-level sufficient statistics) and 'FIML'
  # (wide format), chosen at model creation; the two use different sample
  # statistics, so switching afterwards is not possible (mirrors ml_lvm):
  if (x@model == "ml_varcov" && estimator != x@estimator){
    stop("The estimator of an ml_varcov model cannot be changed after the model is created (the 'ML' and 'FIML' estimators use different sample statistics). Rebuild the model with ml_varcov(..., estimator = '", estimator, "').")
  }

  # Robust ML estimators map internally to estimator = "ML" plus a robust
  # configuration. Resolve the requested name into the internal estimator and
  # the robust config now. MLR additionally supports missing data: when the
  # model carries FIML data (per-observation missingness patterns), the internal
  # estimator is "FIML" rather than "ML" (mirroring the constructor, where
  # estimator = "ML" + missing data auto-switches to "FIML"):
  if (estimator %in% .robustEstimators){
    robust_cfg <- .robust_config_for(estimator)
    if (estimator == "MLR" && length(x@sample@fimldata) > 0){
      internal_estimator <- "FIML"
    } else {
      internal_estimator <- "ML"
    }
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

  # Robust ML requires raw data: the sandwich corrections need the asymptotic
  # covariance Gamma (MLM family) or casewise scores (MLR), neither of which is
  # available from summary statistics alone.
  #   * MLR (Huber-White SEs + Yuan-Bentler-Mplus test) supports BOTH complete
  #     data (estimator = "ML" internally) and within-row missing data (FIML):
  #     the pattern-wise casewise scores and the FIML observed information are
  #     formed from the stored raw data and missingness patterns.
  #   * The MLM family (robust.sem / Satorra-Bentler) requires the asymptotic
  #     covariance Gamma of the COMPLETE-data sample statistics, which is not
  #     defined under missingness; it therefore supports complete data only.
  if (length(robust_cfg) > 0){
    has_raw <- nrow(x@sample@rawdata) > 0 || length(x@sample@fimldata) > 0
    if (!has_raw){
      stop("estimator = '", robust_cfg$label, "' (robust ML) requires raw data, which ",
           "is not stored in this model. Rebuild the model from a raw data frame ",
           "(the `data` argument) before switching to '", robust_cfg$label, "'.")
    }
    has_missing <- length(x@sample@fimldata) > 0
    if (has_missing && robust_cfg$label != "MLR"){
      stop("estimator = '", robust_cfg$label, "' (robust.sem / Satorra-Bentler) ",
           "supports complete data only: it requires the asymptotic covariance ",
           "(Gamma) of the complete-data sample statistics, which is not defined ",
           "under missing data. The model contains missing data. Use estimator = ",
           "'MLR', which provides robust standard errors and a scaled test ",
           "statistic under FIML (missing data).")
    }
    if (has_missing && robust_cfg$label == "MLR" && nrow(x@sample@rawdata) == 0){
      # MLR/FIML needs the stored raw data (with NAs) for the casewise scores.
      stop("estimator = 'MLR' with missing data requires the raw data to be ",
           "stored (storedata = TRUE). Rebuild the model from a raw data frame.")
    }
  }

  # Ensure the asymptotic covariance Gamma of the sample statistics is present,
  # computing it from the stored raw data if missing. Needed by the robust.sem
  # estimators (MLM/MLMV/MLMVS); since 0.16.5 Gamma is no longer computed for
  # plain ML fits (V2-3), so switching such a model to a robust.sem estimator
  # must (re)compute it here -- otherwise the robust sandwich SEs would silently
  # fall back to the naive VCOV.
  ensure_ml_gamma <- function(model){
    if (has_WLS_Gamma(model)) return(model)
    if (!.hasSlot(model@sample, "rawdata")) return(model)
    rd <- model@sample@rawdata
    vars <- attr(rd, "vars"); groupcol <- attr(rd, "groups")
    if (is.null(vars) || is.null(groupcol) || nrow(rd) == 0) return(model)
    nG <- nrow(model@sample@groups); labs <- model@sample@groups$label
    gm <- vector("list", nG); ok <- TRUE
    for (g in seq_len(nG)){
      sub <- as.matrix(rd[rd[[groupcol]] == labs[g], vars, drop = FALSE])
      sub <- sub[stats::complete.cases(sub), , drop = FALSE]
      if (nrow(sub) < 2){ ok <- FALSE; break }
      gm[[g]] <- tryCatch(LS_Gamma(sub, meanstructure = model@meanstructure,
                                   corinput = isTRUE(model@sample@corinput)),
                          error = function(e){ ok <<- FALSE; NULL })
    }
    if (ok) model@sample@WLS.Gamma <- gm
    model
  }

  # Recursively set the estimator on a model and its baseline/saturated:
  .set_one <- function(model, estimator, robust_cfg){
    .check_estimator_requirements(model, estimator)
    model@estimator <- estimator
    if (.hasSlot(model, "robust")){
      model@robust <- robust_cfg
    }
    # Robust ML estimators need Gamma present (robust.sem uses its values;
    # complete-data MLR uses its dimensions). Compute it lazily if a plain-ML
    # model (which since 0.16.5 stores no Gamma, V2-3) is switched to a robust
    # estimator. Skipped for FIML/MLR-missing (no complete-data Gamma):
    if (isTRUE(nzchar(robust_cfg$se)) && estimator != "FIML"){
      model <- ensure_ml_gamma(model)
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
