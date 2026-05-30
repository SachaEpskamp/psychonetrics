# Blume-Capel model creator.
#
# The Blume-Capel model is the full Spin-distribution model: it estimates a
# linear node potential (tau), a quadratic node potential (delta), and pairwise
# interactions (omega):
#
#   P(X = x) propto exp( sum_i tau_i x_i - sum_i delta_i x_i^2 + sum_{i<j} omega_ij x_i x_j )
#
# The Ising model is the special case in which all delta_i are fixed to zero, so
# BlumeCapel() simply dispatches to the shared Ising() engine with the quadratic
# 'delta' parameters freed. By convention the model is defined for the ternary
# -1/0/1 coding (a warning is raised for any other response set), but any
# ordered numeric response set is allowed.
BlumeCapel <- function(
  data, # Dataset
  omega = "full", # Network structure
  tau, # Linear node potentials
  delta, # Quadratic node potentials (default: all free)
  beta,
  beta_model = c("beta","log_beta"),
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  covtype = c("choose","ML","UB"),
  responses, # May not be missing if data is missing. Defaults to -1/0/1.
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer,
  storedata = FALSE,
  WLS.W,
  sampleStats, # Leave to missing
  identify = TRUE,
  verbose = FALSE,
  maxNodes = 20,
  maxStates = 2^maxNodes, # Max number of response patterns enumerated by ML (length(responses)^nNode)
  min_sum = -Inf, # Used for thresholded estimation
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
){
  # Forward every supplied argument verbatim to the shared spinModel() engine,
  # flagging the model type as BlumeCapel (which frees the delta parameters).
  # Using match.call() means arguments the user did not supply stay missing, so
  # the engine's own defaults and missing()-handling (in particular the default
  # free delta) apply unchanged.
  # Use the function object itself (not the name) in the call head: spinModel is
  # an internal, non-exported function and would not be found by name lookup
  # from the caller's environment when the package is installed.
  mc <- match.call()
  mc[[1L]] <- spinModel
  mc$model_name <- "BlumeCapel"
  eval(mc, parent.frame())
}
