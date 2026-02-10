# Penalized ML fit function for Ising distribution
# Wraps the standard ML fit and adds L1 + L2 penalty terms
penMaxLikEstimator_Ising <- function(prep, x, model) {
  if (model@cpp) {
    fit <- maxLikEstimator_Ising_cpp(prep)
  } else {
    fit <- maxLikEstimator_Ising(prep)
  }
  addPenaltyFit(fit, x, model)
}
