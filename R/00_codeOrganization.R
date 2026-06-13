# Numeric code are internal functions:
# 00: This file
# 01: classes definitions
# 03: Functions to aide in model formation. E.g., generate sample table, parameter table, model matrices..
# 04: General Gaussian fit, gradient and Hessian
# 05: General Binary fit, gradient and Hessian
# 06: Model specific fit functions
# 99: Old codes to be removed

# Letter codes are exported functions and executables both used internally and externally (e.g., addfit())
# a: model generation codes
# b: functions that expand models, e.g., add MIs
# c: model execution
# d: Stepup search

# I copied this piece of code from Lavaan mainly:

.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version),"! For questions, issues, and bug reports, please see github.com/SachaEpskamp/psychonetrics.")
}

# Helper to emit a one-time experimental feature warning (only when version < 0.17):
.experimental_warned <- new.env(parent = emptyenv())

# Historic helper (0.15.31 development): the two-level sufficient-statistics
# ML estimator (ml_lvm with estimator = "ML", distribution "TwoLevelGaussian")
# initially only had an R implementation, and this helper downgraded such
# models to the R code path. Since the C++ twins were added (fit, gradient,
# expected Hessian, model Jacobian and prepare), the C++ path is fully
# supported and DEFAULT for these models, exactly as for the other model
# families, so this is now the identity. It is kept (and still called at the
# central dispatch points) so that the R/C++ switch remains easy to force
# centrally if ever needed again:
force_R_path_if_needed <- function(model){
  model
}

# Safe reader for the two-level sufficient statistics: psychonetrics objects
# saved before 0.15.31 lack the 'twolevel' slot in their samplestats, so guard
# with .hasSlot before accessing:
get_twolevel_stats <- function(samplestats){
  if (methods::.hasSlot(samplestats, "twolevel")){
    samplestats@twolevel
  } else {
    list()
  }
}

# TRUE if a model uses the two-level ML estimator (distribution
# "TwoLevelGaussian") AND its data carry within-cluster missingness (so the
# per-pattern missing-data likelihood is used). Such models use the R fit/
# gradient path and numeric Fisher information for SEs (Phase 4):
twolevel_model_has_missing <- function(model){
  if (model@distribution != "TwoLevelGaussian") return(FALSE)
  ts <- get_twolevel_stats(model@sample)
  length(ts) > 0 && any(vapply(ts, function(s) isTRUE(s$missing), logical(1)))
}

experimentalWarning <- function(feature) {
  # Only warn for pre-0.17 versions (the new 0.16 methods are still flagged
  # experimental; raise this once they are considered stable):
  ver <- utils::packageVersion("psychonetrics")
  if (ver >= "0.17") return(invisible())
  # Only warn once per feature per session:
  if (isTRUE(.experimental_warned[[feature]])) return(invisible())
  .experimental_warned[[feature]] <- TRUE
  message("Note: '", feature, "' is experimental in psychonetrics ", ver,
          ". Please report any unexpected behavior to https://github.com/SachaEpskamp/psychonetrics/issues")
}

# Robust ML configuration accessor. The @robust slot was added in 0.15.31;
# models saved before that lack the slot, so it is always read via .hasSlot()
# and treated as a non-robust (empty) configuration if absent. Returns a list
# with components:
#   se    : "robust.sem" (MLM family), "robust.huber.white" (MLR), or "" (none)
#   test  : "satorra.bentler", "scaled.shifted", "mean.var.adjusted",
#           "yuan.bentler.mplus", or "" (none)
#   label : the user-facing estimator name ("MLM"/"MLMV"/"MLMVS"/"MLR")
# Only meaningful when x@estimator == "ML" (robust ML maps internally to ML
# point estimation plus robust SE/test flags).
get_robust_config <- function(x){
  if (!.hasSlot(x, "robust")) return(list())
  cfg <- x@robust
  if (!is.list(cfg)) return(list())
  cfg
}

# TRUE if the model carries a robust ML configuration:
is_robust_ML <- function(x){
  cfg <- get_robust_config(x)
  isTRUE(nzchar(cfg$se)) || isTRUE(nzchar(cfg$test))
}

# Gaussian likelihood scaling for the varcov / lvm families: "normal" (the
# default, n denominator: chisq = N * Fhat, naive VCOV = Info^-1 / N) or
# "wishart" (n-1 denominator: unbiased sample covariance, chisq = (N-1) * Fhat,
# VCOV = Info^-1 / (N-1)), matching lavaan's likelihood = argument. Stored in
# x@types$likelihood by the constructors; read here with a default of "normal"
# so that models created before this feature (no likelihood entry) behave
# exactly as before.
model_likelihood <- function(x){
  lik <- x@types$likelihood
  if (is.null(lik) || !nzchar(lik)) return("normal")
  lik
}

# TRUE if the model uses the Wishart Gaussian likelihood scaling.
is_wishart <- function(x){
  identical(model_likelihood(x), "wishart")
}
