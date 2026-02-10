# Penalty helper functions for Penalized Maximum Likelihood (PML) estimation
# These add L1 (LASSO) and/or L2 (ridge) penalty terms to fit values and gradients

# Add penalty to gradient vector (L2 gradient + L1 subgradient)
addPenaltyGradient <- function(grad, x, model) {
  lambda_vec <- penaltyVector(model)
  alpha <- model@penalty$alpha
  if (is.null(alpha)) alpha <- 1
  # L2 gradient: lambda * (1-alpha) * x
  grad <- grad + lambda_vec * (1 - alpha) * x
  # L1 subgradient: lambda * alpha * sign(x) for x != 0, 0 for x == 0
  grad <- grad + lambda_vec * alpha * sign(x)
  grad
}

# Add penalty to fit value (L2 + L1)
addPenaltyFit <- function(fit, x, model) {
  lambda_vec <- penaltyVector(model)
  alpha <- model@penalty$alpha
  if (is.null(alpha)) alpha <- 1
  # L2 (ridge): 0.5 * sum(lambda * (1-alpha) * x^2)
  # L1 (LASSO): sum(lambda * alpha * |x|)
  fit + 0.5 * sum(lambda_vec * (1 - alpha) * x^2) +
    sum(lambda_vec * alpha * abs(x))
}
