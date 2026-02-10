# Penalized FIML fit function for Gaussian distribution
# Wraps the standard FIML fit and adds L1 + L2 penalty terms
penFIMLEstimator_Gauss <- function(prep, x, model) {
  if (model@cpp) {
    fit <- fimlestimator_Gauss_cpp(prep)
  } else {
    fit <- fimlEstimator_Gauss(prep)
  }
  addPenaltyFit(fit, x, model)
}
