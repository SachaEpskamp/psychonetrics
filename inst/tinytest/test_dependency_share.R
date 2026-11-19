# Tests for dependency_share(): the share of the dependence between composites
# that is carried by a composite model (cvm / cnm).
#
# Most tests use a covariance matrix as input. With one composite of several
# indicators and further only composites with a single indicator, the model is
# the canonical correlation model with one canonical correlation, so that the
# chi-squares of the model and of the baseline model follow from the canonical
# correlations between the two sets of variables. This gives expected values
# that do not depend on the package.

suppressMessages(library(psychonetrics))

quiet <- function(expr) suppressMessages(expr)
share <- function(x) quiet(dependency_share(x, verbose = FALSE))

## ---- Covariance matrix of five variables ----
# y1 and y2 form a composite; y3, y4 and y5 are composites with one indicator.
L <- matrix(c(1.0,  0.0,  0.0, 0.0, 0.0,
              0.5,  0.9,  0.0, 0.0, 0.0,
              0.4,  0.3,  1.1, 0.0, 0.0,
              0.3, -0.2,  0.4, 0.8, 0.0,
              0.2,  0.4, -0.1, 0.3, 1.2), 5, 5, byrow = TRUE)
S <- L %*% t(L)
vars <- paste0("y", 1:5)
dimnames(S) <- list(vars, vars)
N <- 400
A <- 1:2
B <- 3:5
W <- matrix(0, 5, 4)
W[cbind(1:5, c(1, 1, 2, 3, 4))] <- 1
comps <- c("C", "y3", "y4", "y5")

# Expected values from the canonical correlations:
r2 <- sort(Re(eigen(solve(S[A, A], S[A, B]) %*% solve(S[B, B], S[B, A]))$values), decreasing = TRUE)
chisq0 <- -N * sum(log(1 - r2))
chisq <- -N * log(1 - r2[2])
expected <- 1 - (chisq - 2) / (chisq0 - 6)

## ---- 1. Free weights: the model with one canonical correlation ----
mod <- quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars, latents = comps)))
expect_equal(mod@fitmeasures$df, 2)
expect_equal(mod@fitmeasures$chisq, chisq, tolerance = 1e-4)
expect_equal(share(mod), expected, tolerance = 1e-5)

# The sets (composites with a single indicator form one set) and the baseline model:
set <- psychonetrics:::dependency_share_sets(mod)
expect_equal(unname(set), c(1, 1, 0, 0, 0))
base <- quiet(psychonetrics:::dependency_share_baseline(mod, set))
expect_equal(base@fitmeasures$df, 6)
expect_equal(base@fitmeasures$chisq, chisq0, tolerance = 1e-4)

# The share is the RNI of the model when this baseline model replaces the default one:
mod_rni <- mod
mod_rni@baseline_saturated$baseline <- base
mod_rni <- quiet(addfit(mod_rni))
expect_equal(mod_rni@fitmeasures$rni, share(mod), tolerance = 1e-8)

# Printed output and returned value:
out <- capture.output(res <- quiet(dependency_share(mod)))
expect_equal(out, paste0("The calculated share of the dependency equals: ", formatC(expected, digits = 3, format = "f")))
expect_equal(capture.output(res2 <- quiet(dependency_share(mod, verbose = FALSE))), character(0))
expect_false(withVisible(quiet(dependency_share(mod, verbose = FALSE)))$visible)
expect_true(is.numeric(res) && length(res) == 1 && is.null(attributes(res)))
expect_equal(res, res2)

## ---- 2. With a mean structure the share is the same ----
mod_m <- quiet(runmodel(cvm(covs = S, means = c(1, 2, 3, 4, 5), nobs = N, covtype = "ML", weights = W, vars = vars,
                            latents = comps)))
expect_true(mod_m@meanstructure)
expect_equal(share(mod_m), expected, tolerance = 1e-5)

## ---- 3. Fixed weights (sum score) ----
# The other variables then only depend on y1 and y2 through their sum:
cSum <- S[A, A] %*% c(1, 1)
SB_c <- S[B, B] - S[B, A] %*% c(1, 1) %*% t(c(1, 1)) %*% S[A, B] / sum(S[A, A])
SB_A <- S[B, B] - S[B, A] %*% solve(S[A, A], S[A, B])
chisq_sum <- N * (log(det(SB_c)) - log(det(SB_A)))
mod_sum <- quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars, latents = comps,
                              sumscores = TRUE)))
expect_equal(mod_sum@fitmeasures$df, 3)
expect_equal(mod_sum@fitmeasures$chisq, chisq_sum, tolerance = 1e-4)
expect_equal(share(mod_sum), 1 - (chisq_sum - 3) / (chisq0 - 6), tolerance = 1e-5)

## ---- 4. Other types of the composites and of their relations ----
types <- if (at_home()) list(c("ggm", "ggm"), c("cov", "ggm"), c("chol", "cov"), c("prec", "cov"), c("cor", "cov")) else
  list(c("ggm", "ggm"))
for (type in types){
  mod_t <- quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars, latents = comps,
                              composite = type[1], latent = type[2])))
  expect_equal(share(mod_t), expected, tolerance = 1e-5, info = paste(type, collapse = " / "))
}

## ---- 5. Two composites with several indicators and two with one ----
L6 <- matrix(c(1.0,  0.0,  0.0,  0.0, 0.0, 0.0,
               0.5,  0.9,  0.0,  0.0, 0.0, 0.0,
               0.4,  0.3,  1.1,  0.0, 0.0, 0.0,
               0.3, -0.2,  0.4,  0.8, 0.0, 0.0,
               0.2,  0.4, -0.1,  0.3, 1.2, 0.0,
               0.3,  0.1,  0.2, -0.2, 0.4, 1.0), 6, 6, byrow = TRUE)
S6 <- L6 %*% t(L6)
vars6 <- paste0("y", 1:6)
dimnames(S6) <- list(vars6, vars6)
W6 <- matrix(0, 6, 4)
W6[cbind(1:6, c(1, 1, 2, 2, 3, 4))] <- 1
mod6 <- quiet(runmodel(cnm(covs = S6, nobs = N, covtype = "ML", weights = W6, vars = vars6)))
set6 <- psychonetrics:::dependency_share_sets(mod6)
expect_equal(unname(set6), c(1, 1, 2, 2, 0, 0))
base6 <- quiet(psychonetrics:::dependency_share_baseline(mod6, set6))
# Covariances between the three sets: 15 - 3 within the sets
expect_equal(base6@fitmeasures$df, 12)
# Closed form of the chi-square of a model with independent sets:
chisq0_6 <- N * (log(det(S6[1:2, 1:2])) + log(det(S6[3:4, 3:4])) + log(det(S6[5:6, 5:6])) - log(det(S6)))
expect_equal(base6@fitmeasures$chisq, chisq0_6, tolerance = 1e-4)
expect_equal(share(mod6), 1 - (mod6@fitmeasures$chisq - mod6@fitmeasures$df) / (chisq0_6 - 12), tolerance = 1e-5)

## ---- 6. Models for which the share is not available ----
expect_error(dependency_share(1), "cvm\\(\\) or cnm\\(\\)")
vc <- quiet(runmodel(varcov(covs = S, nobs = N, covtype = "ML", vars = vars)))
expect_error(share(vc), "only available for composite models")
expect_error(share(quiet(cvm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars))), "runmodel")

# Fixed edge or covariance between the composites, as after pruning:
mod_cnm <- quiet(cnm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars, latents = comps))
mod_fix <- quiet(runmodel(fixpar(mod_cnm, "omega_zeta", row = 3, col = 2)))
expect_error(share(mod_fix), "1 of the 6 edges between the composites \\('omega_zeta'\\) are fixed")
mod_fix2 <- quiet(runmodel(fixpar(mod, "sigma_zeta", row = 4, col = 3)))
expect_error(share(mod_fix2), "1 of the 6 covariances between the composites \\('sigma_zeta'\\) are fixed")

# Edges constrained to be equal:
edges <- which(mod_cnm@parameters$matrix == "omega_zeta" & !mod_cnm@parameters$fixed)
mod_eq <- quiet(runmodel(parequal(mod_cnm, edges[1:2])))
expect_equal(mod_eq@fitmeasures$df, 3)
expect_error(share(mod_eq), "constrained to be equal")

# Fixed element within a composite:
mod_within <- quiet(runmodel(fixpar(mod, "sigma_c", row = 2, col = 1)))
expect_error(share(mod_within), "within the composites is not saturated: 1 of the 6 elements")

# Effect between the composites:
beta <- matrix(0, 4, 4)
beta[2, 1] <- 1
mod_beta <- quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = W, beta = beta, vars = vars,
                               latents = comps)))
expect_error(share(mod_beta), "'beta'")

# Common factor:
lambda <- matrix(0, 6, 3)
lambda[3:4, 2] <- 1
W6f <- matrix(0, 6, 3)
W6f[1:2, 1] <- 1
W6f[5:6, 3] <- 1
mod_fac <- suppressWarnings(quiet(runmodel(cvm(covs = S6, nobs = N, covtype = "ML", weights = W6f, lambda = lambda,
                                               vars = vars6, latents = c("C1", "F", "C2")))))
expect_error(share(mod_fac), "1 common factor\\(s\\) \\(F\\)")

# Other estimator:
mod_uls <- suppressWarnings(quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = W, vars = vars,
                                               latents = comps, estimator = "ULS"))))
expect_error(share(mod_uls), "only available for the 'ML' and 'FIML' estimators")

# Multiple groups:
mod_groups <- quiet(runmodel(cvm(covs = list(S, S), nobs = c(N, N), covtype = "ML", weights = W, vars = vars,
                                 latents = comps)))
expect_error(share(mod_groups), "multiple groups are not supported")

# Only composites with a single indicator (a saturated model), or a single composite:
mod_single <- quiet(runmodel(cvm(covs = S, nobs = N, covtype = "ML", weights = diag(5), vars = vars)))
expect_error(share(mod_single), "all composites have a single indicator")

# A different saturated model than the baseline model:
mod_sat <- mod
mod_sat@fitmeasures$unrestricted.logl <- mod@fitmeasures$unrestricted.logl + 1
expect_error(share(mod_sat), "saturated log-likelihood")

# Composites that are (nearly) unrelated:
S0 <- S
S0[A, B] <- 0.01 * S[A, B]
S0[B, A] <- t(S0[A, B])
mod_0 <- suppressWarnings(quiet(runmodel(cvm(covs = S0, nobs = N, covtype = "ML", weights = W, vars = vars,
                                             latents = comps))))
expect_error(share(mod_0), "not clearly related")

## ---- 7. Raw data with missing values (FIML) ----
data("NA2020")
vars_na <- c("regular_sleep", "worried_sleep", "sleep_interfere", "happy_health", "optimistic_future", "very_happy",
             "feel_alone", "happy_love_life")
W_na <- matrix(0, 8, 7)
W_na[cbind(1:8, c(1, 2, 2, 3, 4, 5, 6, 7))] <- 1
mod_na <- quiet(runmodel(cnm(NA2020, weights = W_na, vars = vars_na, estimator = "FIML")))
set_na <- psychonetrics:::dependency_share_sets(mod_na)
expect_equal(unname(set_na), c(0, 2, 2, 0, 0, 0, 0, 0))
base_na <- quiet(psychonetrics:::dependency_share_baseline(mod_na, set_na))
# The two sleep problems against the six other variables:
expect_equal(base_na@fitmeasures$df, 12)
expect_equal(base_na@fitmeasures$unrestricted.logl, mod_na@fitmeasures$unrestricted.logl, tolerance = 1e-8)
expect_equal(base_na@fitmeasures$chisq, 89.4632, tolerance = 1e-4)
expect_equal(mod_na@fitmeasures$df, 5)
expect_equal(share(mod_na), 1.015036, tolerance = 1e-5)
# Both weights of the composite fixed to 1:
mod_na1 <- quiet(runmodel(cnm(NA2020, weights = W_na, vars = vars_na, estimator = "FIML", sumscores = TRUE)))
expect_equal(mod_na1@fitmeasures$df, 6)
expect_equal(share(mod_na1), 1.016229, tolerance = 1e-5)
