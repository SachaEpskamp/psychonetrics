# Proximal Gradient (ISTA) optimizer for Penalized ML estimation
# Handles L1 (LASSO) and elastic net penalties via soft-thresholding

psychonetrics_proximal_gradient <- function(x, lower, upper, bounded, tol = 1e-8) {
  # Extract parameters
  par <- parVector(x)
  nPar <- length(par)

  if (nPar == 0) {
    x@computed <- TRUE
    x@optim <- list(par = par, value = 0, convergence = 0,
                    message = "No free parameters", optimizer = "proximal_gradient")
    return(x)
  }

  lambda_vec <- penaltyVector(x)  # per-free-parameter penalty lambda
  alpha <- x@penalty$alpha
  if (is.null(alpha)) alpha <- 1  # default to LASSO

  verbose <- x@verbose

  # ISTA settings
  max_iter <- 10000
  # tol is now a function parameter (default 1e-8)
  t_init <- 1.0     # initial step size
  beta_bt <- 0.5    # backtracking reduction factor
  t_min <- 1e-14    # minimum step size
  fncount <- 0
  grcount <- 0

  # Soft-thresholding operator
  soft_thresh <- function(z, thresh) {
    sign(z) * pmax(abs(z) - thresh, 0)
  }

  # Choose fit/gradient functions based on cpp flag
  if (x@cpp) {
    fn <- function(par) psychonetrics_fitfunction_cpp(par, x)
    gr <- function(par) psychonetrics_gradient_cpp(par, x)
  } else {
    fn <- function(par) psychonetrics_fitfunction(par, x)
    gr <- function(par) psychonetrics_gradient(par, x)
  }

  # NOTE: fn() now returns ML + L2 + L1 (full penalized objective)
  # NOTE: gr() now returns ML_grad + L2_grad + L1_subgrad (full subgradient)
  # ISTA needs smooth vs non-smooth separation, so we subtract the L1 parts.

  # Smooth objective = fn(par) - L1 = (ML + L2 + L1) - L1 = ML + L2
  smooth_obj <- function(par) {
    full_fit <- fn(par)
    fncount <<- fncount + 1
    if (!is.finite(full_fit)) return(full_fit)
    l1 <- sum(lambda_vec * alpha * abs(par))
    full_fit - l1
  }

  # Full penalized objective (for convergence checking)
  full_obj <- function(par, smooth_val) {
    l1 <- sum(lambda_vec * alpha * abs(par))
    smooth_val + l1
  }

  # Smooth gradient = gr(par) - L1_subgrad = (ML_grad + L2_grad + L1_subgrad) - L1_subgrad
  smooth_grad <- function(par) {
    full_grad <- as.vector(gr(par))
    grcount <<- grcount + 1
    # Subtract L1 subgradient to get smooth part only
    l1_subgrad <- lambda_vec * alpha * sign(par)
    full_grad - l1_subgrad
  }

  # Compute initial objective
  f_x <- smooth_obj(par)
  if (!is.finite(f_x)) {
    stop("Initial parameter values produce non-finite objective. Try different starting values.")
  }
  obj_x <- full_obj(par, f_x)

  convergence <- 0
  message_out <- "converged"
  t_step <- t_init

  for (iter in seq_len(max_iter)) {
    # Compute gradient of smooth part
    grad <- smooth_grad(par)

    # Backtracking line search
    t_k <- t_step
    step_found <- FALSE

    while (t_k >= t_min) {
      # Gradient step
      z <- par - t_k * grad

      # Proximal step: soft-threshold on penalized parameters
      par_new <- z
      pen_idx <- which(lambda_vec > 0)
      if (length(pen_idx) > 0) {
        par_new[pen_idx] <- soft_thresh(z[pen_idx], t_k * lambda_vec[pen_idx] * alpha)
        # Elastic net scaling (denominator adjustment for L2)
        if (alpha < 1) {
          par_new[pen_idx] <- par_new[pen_idx] / (1 + t_k * lambda_vec[pen_idx] * (1 - alpha))
        }
      }

      # Enforce bounds
      if (bounded) {
        par_new <- pmax(pmin(par_new, upper), lower)
      }

      # Evaluate smooth part at proposed point
      f_new <- smooth_obj(par_new)

      if (!is.finite(f_new)) {
        t_k <- t_k * beta_bt
        next
      }

      # Sufficient decrease condition (proximal Armijo):
      # f_smooth(x_new) <= f_smooth(x) + grad' * (x_new - x) + (1/(2*t_k)) * ||x_new - x||^2
      diff_vec <- par_new - par
      quad_approx <- f_x + sum(grad * diff_vec) + sum(diff_vec^2) / (2 * t_k)

      if (f_new <= quad_approx + 1e-12) {
        step_found <- TRUE
        f_x <- f_new
        break
      }

      t_k <- t_k * beta_bt
    }

    if (!step_found) {
      convergence <- 1
      message_out <- "Backtracking failed: step size too small"
      if (verbose) message(paste0("ISTA iter ", iter, ": backtracking failed"))
      break
    }

    # Compute new full objective
    obj_new <- full_obj(par_new, f_x)

    # Check convergence on relative objective change
    rel_change <- abs(obj_x - obj_new) / (abs(obj_x) + 1e-10)

    # Update
    par <- par_new
    obj_x <- obj_new

    # Try increasing step size for next iteration (adaptive restart)
    t_step <- min(t_k / beta_bt, t_init)

    if (rel_change < tol) {
      if (verbose) message(paste0("ISTA converged at iteration ", iter,
                                   " (rel_change = ", formatC(rel_change, format = "e", digits = 2), ")"))
      break
    }

    # Periodic progress messages
    if (verbose && iter %% 100 == 0) {
      n_zero <- sum(abs(par[lambda_vec > 0]) < 1e-10)
      n_pen <- sum(lambda_vec > 0)
      message(paste0("ISTA iter ", iter, ": obj = ", formatC(obj_x, digits = 6),
                      ", zeros = ", n_zero, "/", n_pen,
                      ", step = ", formatC(t_k, format = "e", digits = 2)))
    }
  }

  if (iter >= max_iter) {
    convergence <- 1
    message_out <- "Maximum iterations reached"
    if (verbose) message("ISTA: maximum iterations reached")
  }

  # Summary stats
  n_zero <- sum(abs(par[lambda_vec > 0]) < 1e-10)
  n_pen <- sum(lambda_vec > 0)

  # Update model and return
  x <- updateModel(par, x, updateMatrices = TRUE)
  x@computed <- TRUE
  x@objective <- obj_x
  x@optim <- list(
    par = par,
    value = obj_x,
    convergence = convergence,
    message = message_out,
    fncount = fncount,
    grcount = grcount,
    optimizer = "proximal_gradient",
    penalty = list(
      n_penalized = n_pen,
      n_zero = n_zero,
      lambda = x@penalty$lambda,
      alpha = alpha
    )
  )

  x
}
