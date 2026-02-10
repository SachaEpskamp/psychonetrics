# Post-selection refit: fix penalized zeros, remove penalty, re-run with standard ML
# This allows obtaining standard errors and fit measures after penalized estimation

refit <- function(
    x,                 # A computed penalized psychonetrics model
    threshold = 1e-8,  # Threshold below which penalized parameters are considered zero
    verbose,
    log = TRUE,
    ...
) {
  if (!is(x, "psychonetrics")) stop("input must be a psychonetrics object")
  if (!x@computed) stop("Model must be computed before refitting. Run runmodel() first.")
  if (missing(verbose)) verbose <- x@verbose

  # Check if model has penalty
  pen_lambda <- x@penalty$lambda
  if (is.null(pen_lambda) || pen_lambda == 0) {
    if (all(x@parameters$penalty_lambda == 0)) {
      message("Model has no penalty. Nothing to refit.")
      return(x)
    }
  }

  # Identify penalized parameters that are (near) zero
  parvec <- parVector(x)
  lambdas <- penaltyVector(x)
  zero_pars <- which(lambdas > 0 & abs(parvec) < threshold)

  if (verbose) {
    message(paste0("Fixing ", length(zero_pars), " of ", sum(lambdas > 0),
                    " penalized parameters to zero for refit"))
  }

  # Create new model: fix the zero parameters
  new_mod <- x

  if (length(zero_pars) > 0) {
    for (par_idx in zero_pars) {
      rows <- which(new_mod@parameters$par == par_idx)
      new_mod@parameters$est[rows] <- 0
      new_mod@parameters$par[rows] <- 0
      new_mod@parameters$fixed[rows] <- TRUE
    }

    # Clear derived columns for fixed parameters
    new_mod@parameters <- clearpars(new_mod@parameters,
                                     which(new_mod@parameters$par == 0 &
                                             new_mod@parameters$penalty_lambda > 0))

    # Relabel free parameter indices
    new_mod@parameters <- parRelabel(new_mod@parameters)
  }

  # Remove penalty: set all penalty_lambda to 0
  new_mod@parameters$penalty_lambda <- 0
  new_mod@penalty <- list(lambda = 0, alpha = 1)

  # Switch estimator back to ML (or FIML for FIPML)
  new_mod@estimator <- if (x@estimator == "FIPML") "FIML" else "ML"
  new_mod@optimizer <- "default"
  new_mod@computed <- FALSE

  # Log entry
  if (log) {
    new_mod <- addLog(new_mod, paste0("Post-selection refit: fixed ",
                                       length(zero_pars), " penalized parameters to zero"))
  }

  # Re-run the model with standard ML
  new_mod <- runmodel(new_mod, verbose = verbose, ...)

  new_mod
}
