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
  packageStartupMessage("This is ",paste(pkgname, version),"! Note: this is BETA software! Please mind that the package may not be stable and report any bugs! For more information, please see psychonetrics.org, for questions and issues, please see github.com/SachaEpskamp/psychonetrics.")
}
