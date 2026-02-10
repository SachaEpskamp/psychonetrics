# Automated lambda selection for PML/PFIML estimation via EBIC grid search
# Implements: find_penalized_lambda(), compute_lambda_max(), compute_penalized_ebic(),
#             apply_beta_min(), parse_ebic_gamma()

# Parse a criterion string into an EBIC gamma value
# Returns NULL for non-EBIC criteria (fall back to addfit)
parse_ebic_gamma <- function(criterion) {
  crit <- tolower(criterion)
  if (crit == "bic") return(0)
  if (crit == "ebic.25") return(0.25)
  if (crit == "ebic.5") return(0.5)
  if (crit == "ebic.75") return(0.75)  # Note: addfit uses gamma=0.7 for this label
  if (crit == "ebic1") return(1.0)
  NULL
}


# Compute the analytical lambda_max from KKT conditions
# This is the smallest lambda for which all penalized parameters are exactly zero.
# Derived from: lambda_max = max(|grad_j|) / alpha  at the null-penalized model.
compute_lambda_max <- function(model) {
  alpha <- model@penalty$alpha
  if (is.null(alpha)) alpha <- 1

  # Identify which parameter table rows are penalized (have NA lambda = auto-select)
  is_pen_row <- is.na(model@parameters$penalty_lambda) &
                model@parameters$par > 0 &
                !model@parameters$fixed

  # Create the null model: fix all penalized params to 0, use ML, no penalty
  null_mod <- model
  pen_par_ids <- unique(null_mod@parameters$par[is_pen_row])

  for (pid in pen_par_ids) {
    rows <- which(null_mod@parameters$par == pid)
    null_mod@parameters$est[rows] <- 0
    null_mod@parameters$par[rows] <- 0
    null_mod@parameters$fixed[rows] <- TRUE
  }

  # Remove penalty and switch to ML/FIML
  null_mod@parameters$penalty_lambda <- 0
  null_mod@penalty <- list(lambda = 0, alpha = 1)
  null_mod@estimator <- if (model@estimator == "PFIML") "FIML" else "ML"
  null_mod@optimizer <- "default"
  null_mod@computed <- FALSE

  # Relabel free parameters
  null_mod@parameters <- parRelabel(null_mod@parameters)

  # Fit the null model
  null_mod <- suppressWarnings(suppressMessages(
    runmodel(null_mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE,
             addInformation = FALSE, verbose = FALSE, log = FALSE,
             warn_improper = FALSE, warn_gradient = FALSE, warn_bounds = FALSE)
  ))

  # Now build a temporary full model at the null solution:
  # All params present (as in original model), penalized ones at 0, no penalty
  temp_mod <- model
  temp_mod@parameters$penalty_lambda <- 0
  temp_mod@penalty <- list(lambda = 0, alpha = 1)
  temp_mod@estimator <- if (model@estimator == "PFIML") "FIML" else "ML"

  # Transfer non-penalized parameter estimates from null model
  # Match by (matrix, row, col, group_id) since par indices changed after parRelabel
  for (pid in unique(null_mod@parameters$par[null_mod@parameters$par > 0])) {
    null_rows <- which(null_mod@parameters$par == pid)
    nr <- null_rows[1]
    # Find matching row in temp model
    match_idx <- which(
      temp_mod@parameters$matrix == null_mod@parameters$matrix[nr] &
      temp_mod@parameters$row == null_mod@parameters$row[nr] &
      temp_mod@parameters$col == null_mod@parameters$col[nr] &
      temp_mod@parameters$group_id == null_mod@parameters$group_id[nr]
    )
    if (length(match_idx) > 0) {
      temp_mod@parameters$est[match_idx] <- null_mod@parameters$est[nr]
    }
  }

  # Set penalized parameters to 0 in temp model
  temp_mod@parameters$est[is_pen_row] <- 0

  # Compute gradient at the null solution in the full model
  if (temp_mod@cpp) {
    grad <- psychonetrics_gradient_cpp(parVector(temp_mod), temp_mod)
  } else {
    grad <- psychonetrics_gradient(parVector(temp_mod), temp_mod)
  }

  # Extract gradient elements corresponding to penalized parameters
  pen_par_vec <- penaltyVector(model)  # Returns 0 for NA (auto-select)
  # We need the par indices that are penalized (had NA lambda)
  pen_free_ids <- unique(model@parameters$par[is_pen_row])

  grad_pen <- abs(grad[pen_free_ids])

  lambda_max <- max(grad_pen) / alpha

  lambda_max
}


# Compute EBIC for a penalized model (Convention A: unpenalized log-likelihood)
compute_penalized_ebic <- function(model, gamma = 0.5, beta_min = 0) {
  # Compute unpenalized log-likelihood at current estimates
  LL <- psychonetrics_logLikelihood(model)

  # Count parameters:
  # npar = (unpenalized free params) + (penalized params with |est| >= beta_min)
  par_table <- model@parameters
  is_pen <- !is.na(par_table$penalty_lambda) & par_table$penalty_lambda > 0 & !par_table$fixed
  is_free_unpen <- par_table$par > 0 & !par_table$fixed &
                   (is.na(par_table$penalty_lambda) | par_table$penalty_lambda == 0)

  # Count unique free unpenalized parameters
  n_unpen <- length(unique(par_table$par[is_free_unpen]))

  # Count unique penalized parameters with |est| >= beta_min
  pen_par_ids <- unique(par_table$par[is_pen])
  n_nonzero_pen <- 0
  for (pid in pen_par_ids) {
    est_val <- par_table$est[which(par_table$par == pid)[1]]
    if (abs(est_val) >= beta_min) {
      n_nonzero_pen <- n_nonzero_pen + 1
    }
  }

  npar <- n_unpen + n_nonzero_pen

  # Sample size and number of variables
  n <- sum(model@sample@groups$nobs)
  nVar <- nrow(model@sample@variables)

  # EBIC = -2*LL + npar*log(n) + 4*npar*gamma*log(nVar)
  ebic <- -2 * LL + npar * log(n) + 4 * npar * gamma * log(nVar)

  ebic
}


# Apply beta_min threshold: set penalized parameters with |est| < beta_min to zero
apply_beta_min <- function(model, beta_min) {
  par_table <- model@parameters
  is_pen <- !is.na(par_table$penalty_lambda) & par_table$penalty_lambda > 0 & !par_table$fixed

  pen_par_ids <- unique(par_table$par[is_pen])
  n_thresholded <- 0

  for (pid in pen_par_ids) {
    rows <- which(par_table$par == pid)
    est_val <- par_table$est[rows[1]]
    if (abs(est_val) < beta_min) {
      model@parameters$est[rows] <- 0
      n_thresholded <- n_thresholded + 1
    }
  }

  # Update model matrices to reflect the zeroed estimates
  model <- updateModel(parVector(model), model, updateMatrices = TRUE)

  attr(model, "n_thresholded") <- n_thresholded
  model
}


# Main function: find optimal lambda via grid search
find_penalized_lambda <- function(
    model,
    criterion = "ebic.5",
    nLambda = 50,
    lambda_min_ratio = 1e-4,
    beta_min,            # default: sqrt(log(p)/n)
    patience = 10,       # early stopping: stop after this many non-improving lambdas
    verbose
) {
  if (!is(model, "psychonetrics")) stop("input must be a psychonetrics object")
  if (!model@estimator %in% c("PML", "PFIML")) {
    stop("find_penalized_lambda() requires a PML or PFIML model")
  }
  if (missing(verbose)) verbose <- model@verbose

  # Parse criterion
  gamma <- parse_ebic_gamma(criterion)
  if (is.null(gamma)) {
    stop("Criterion '", criterion, "' is not supported for lambda search. ",
         "Use one of: 'bic', 'ebic.25', 'ebic.5', 'ebic.75', 'ebic1'")
  }

  # Identify penalized parameters (NA lambda = auto-select)
  is_pen_row <- is.na(model@parameters$penalty_lambda) &
                model@parameters$par > 0 &
                !model@parameters$fixed
  pen_par_ids <- unique(model@parameters$par[is_pen_row])
  p <- length(pen_par_ids)
  n <- sum(model@sample@groups$nobs)

  if (p == 0) {
    if (verbose) message("No auto-select (NA) penalty parameters found. Nothing to search.")
    return(model)
  }

  # Default beta_min
  if (missing(beta_min)) {
    beta_min <- sqrt(log(p) / n)
  }

  if (verbose) {
    message(paste0("Lambda search: ", p, " penalized parameters, n = ", n,
                   ", beta_min = ", formatC(beta_min, format = "e", digits = 3),
                   ", criterion = ", criterion))
  }

  # Step 1: Compute lambda_max
  if (verbose) message("Computing lambda_max from KKT conditions...")
  lambda_max <- compute_lambda_max(model)

  if (!is.finite(lambda_max) || lambda_max <= 0) {
    warning("Could not compute a valid lambda_max. Using lambda_max = 1.")
    lambda_max <- 1
  }

  if (verbose) message(paste0("Lambda_max = ", formatC(lambda_max, format = "e", digits = 3)))

  # Step 2: Build log-spaced grid from lambda_max to lambda_min
  lambda_min <- lambda_max * lambda_min_ratio
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_min), length.out = nLambda))

  # Step 3: Bounds for ISTA
  lower <- lowerBound(model)
  upper <- upperBound(model)

  # Step 4: Loop over lambda grid (warm starts)
  best_crit <- Inf
  best_lambda <- lambda_grid[1]
  best_model <- NULL
  current_model <- model
  patience_counter <- 0

  alpha <- model@penalty$alpha
  if (is.null(alpha)) alpha <- 1

  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]

    # Set lambda for all auto-select (NA) rows
    current_model@parameters$penalty_lambda[is_pen_row] <- lam
    current_model@penalty$lambda <- lam
    current_model@computed <- FALSE

    # Update C++ workspace penalty vector without full rebuild
    updateWorkspacePenaltyLambda(penaltyVector(current_model), current_model)

    # Warm start: use previous solution (already in current_model from last iteration)
    # Run ISTA with loose tolerance (re-optimize best model later)
    tryres <- try(suppressWarnings(suppressMessages({
      current_model <- psychonetrics_proximal_gradient(current_model, lower, upper,
                                                       bounded = TRUE, tol = 1e-5)
    })), silent = TRUE)

    if (is(tryres, "try-error")) {
      if (verbose) message(paste0("  Lambda ", i, "/", nLambda,
                                   " (", formatC(lam, format = "e", digits = 3),
                                   "): ISTA failed, skipping"))
      next
    }

    # Compute EBIC
    crit_val <- tryCatch(
      compute_penalized_ebic(current_model, gamma = gamma, beta_min = beta_min),
      error = function(e) Inf
    )

    # Count non-zero penalized params (after beta_min threshold for counting)
    pen_est <- sapply(pen_par_ids, function(pid) {
      current_model@parameters$est[which(current_model@parameters$par == pid)[1]]
    })
    n_nonzero <- sum(abs(pen_est) >= beta_min)
    n_zero <- p - n_nonzero

    if (verbose && (i %% 10 == 1 || i == nLambda)) {
      message(paste0("  Lambda ", i, "/", nLambda,
                     " (", formatC(lam, format = "e", digits = 3),
                     "): ", criterion, " = ", formatC(crit_val, digits = 6),
                     ", edges = ", n_nonzero, "/", p))
    }

    if (is.finite(crit_val) && crit_val < best_crit) {
      best_crit <- crit_val
      best_lambda <- lam
      best_model <- current_model
      patience_counter <- 0
    } else {
      patience_counter <- patience_counter + 1
      if (patience_counter >= patience) {
        if (verbose) message(paste0("  Early stopping at lambda ", i, "/", nLambda,
                                     ": no improvement for ", patience, " consecutive lambdas"))
        break
      }
    }
  }

  if (is.null(best_model)) {
    warning("Lambda search failed: no valid model found across the grid. Returning original model.")
    return(model)
  }

  # Step 5: Re-optimize best model with tight tolerance
  if (verbose) message("Re-optimizing best model with tight convergence tolerance...")
  best_model@computed <- FALSE
  updateWorkspacePenaltyLambda(penaltyVector(best_model), best_model)
  best_model <- psychonetrics_proximal_gradient(best_model, lower, upper,
                                                 bounded = TRUE, tol = 1e-8)
  # Recompute EBIC with tight solution
  best_crit <- tryCatch(
    compute_penalized_ebic(best_model, gamma = gamma, beta_min = beta_min),
    error = function(e) best_crit
  )

  # Step 6: Apply beta_min threshold on the best model (actual zeroing of estimates)
  best_model <- apply_beta_min(best_model, beta_min)
  n_thresholded <- attr(best_model, "n_thresholded")
  attr(best_model, "n_thresholded") <- NULL

  # Step 6: Store search metadata
  best_model@optim$lambda_search <- list(
    criterion = criterion,
    gamma = gamma,
    best_lambda = best_lambda,
    best_criterion = best_crit,
    beta_min = beta_min,
    lambda_max = lambda_max,
    nLambda = nLambda,
    n_thresholded = n_thresholded
  )

  # Update global penalty lambda
  best_model@penalty$lambda <- best_lambda

  if (verbose) {
    n_pen_final <- sum(!is.na(best_model@parameters$penalty_lambda) &
                       best_model@parameters$penalty_lambda > 0 &
                       !best_model@parameters$fixed &
                       abs(best_model@parameters$est) >= beta_min)
    # Count unique par ids
    pen_rows <- which(!is.na(best_model@parameters$penalty_lambda) &
                      best_model@parameters$penalty_lambda > 0 &
                      !best_model@parameters$fixed)
    pen_pars_final <- unique(best_model@parameters$par[pen_rows])
    n_nonzero_final <- sum(sapply(pen_pars_final, function(pid) {
      abs(best_model@parameters$est[which(best_model@parameters$par == pid)[1]]) >= beta_min
    }))

    message(paste0("\nLambda search complete:",
                   "\n  Best lambda: ", formatC(best_lambda, format = "e", digits = 3),
                   "\n  ", criterion, ": ", formatC(best_crit, digits = 6),
                   "\n  Non-zero penalized: ", n_nonzero_final, "/", p,
                   "\n  Beta_min thresholded: ", n_thresholded, " parameters"))
  }

  best_model
}
