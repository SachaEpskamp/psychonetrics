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

# Helper to emit a one-time experimental feature warning (only when version < 0.15):
.experimental_warned <- new.env(parent = emptyenv())

experimentalWarning <- function(feature) {
  # Only warn for pre-0.15 versions:
  ver <- utils::packageVersion("psychonetrics")
  if (ver >= "0.15") return(invisible())
  # Only warn once per feature per session:
  if (isTRUE(.experimental_warned[[feature]])) return(invisible())
  .experimental_warned[[feature]] <- TRUE
  message("Note: '", feature, "' is experimental in psychonetrics ", ver,
          ". Please report any unexpected behavior to https://github.com/SachaEpskamp/psychonetrics/issues")
}
