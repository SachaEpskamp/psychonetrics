# Penalized ML fit function for Gaussian distribution
# Wraps the standard ML fit and adds L1 + L2 penalty terms
penMaxLikEstimator_Gauss <- function(prep, x, model) {
  if (model@cpp) {
    fit <- maxLikEstimator_Gauss_cpp(prep)
  } else {
    fit <- maxLikEstimator_Gauss(prep)
  }
  addPenaltyFit(fit, x, model)
}
