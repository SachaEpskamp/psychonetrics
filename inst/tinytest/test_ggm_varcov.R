# Tests derived from psychonetrics.org examples 1-4 and 7 (varcov/ggm family),
# reduced to 10 bfi items for speed. Full 25-item versions run at_home() only.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12). Tolerances absorb
# optimizer/BLAS noise across platforms.

suppressMessages(library(psychonetrics))
if (!requireNamespace("psych", quietly = TRUE)) exit_file("psych not available")

bfi <- psych::bfi

## ---- Example 1: saturated GGM (reduced: items A1-A5, C1-C5) ----
bfi10 <- na.omit(bfi[, 1:10])
expect_equal(nrow(bfi10), 2632)

mod_saturated <- suppressWarnings(runmodel(ggm(bfi10)))
expect_inherits(mod_saturated, "psychonetrics")
expect_true(mod_saturated@computed)
expect_true(abs(mod_saturated@fitmeasures$logl - (-42385.3817)) < 0.01)
expect_equal(mod_saturated@fitmeasures$npar, 65)
expect_equal(mod_saturated@fitmeasures$df, 0)

omega <- getmatrix(mod_saturated, "omega")
expect_equal(dim(omega), c(10L, 10L))
expect_true(isSymmetric(unname(omega)))
expect_true(all(diag(omega) == 0))

# Thresholded network has no more edges than the full one:
omega_thr <- getmatrix(mod_saturated, "omega", threshold = TRUE, alpha = 0.01)
expect_equal(dim(omega_thr), c(10L, 10L))
expect_true(sum(omega_thr != 0) <= sum(omega != 0))

# fit() and parameters() invisibly return data frames (printed output captured):
invisible(capture.output({
  ft <- fit(mod_saturated)
  pars <- parameters(mod_saturated)
}))
expect_true(is.data.frame(ft))
expect_true(is.data.frame(pars))

## ---- Example 2: prune + stepup model search ----
mod_pruned <- suppressWarnings(prune(mod_saturated, alpha = 0.01, adjust = "fdr"))
expect_equal(mod_pruned@fitmeasures$df, 23)
expect_true(abs(mod_pruned@fitmeasures$chisq - 72.52659) < 0.01)

mod_stepup <- suppressWarnings(stepup(mod_pruned))
expect_equal(mod_stepup@fitmeasures$df, 20)

cmp <- compare(saturated = mod_saturated, pruned = mod_pruned, stepup = mod_stepup)
expect_true(is.data.frame(cmp))
expect_equal(nrow(cmp), 3)
expect_true(all(is.finite(cmp$AIC)))
expect_true(all(is.finite(cmp$BIC)))
expect_true(all(is.finite(cmp$DF)))
# stepup (df 20) is nested between saturated (df 0) and pruned (df 23):
expect_equal(sort(cmp$DF), c(0, 20, 23))
# The sparser stepup model improves AIC over the saturated model:
expect_true(cmp$AIC[cmp$DF == 20] < cmp$AIC[cmp$DF == 0])

## ---- Example 3: train/test confirmatory GGM ----
set.seed(123)
indices <- sample(nrow(bfi10))
train_data <- bfi10[indices[1:1000], ]
test_data  <- bfi10[indices[1001:nrow(bfi10)], ]
mod_train <- suppressWarnings(
  stepup(prune(runmodel(ggm(train_data)), alpha = 0.01, adjust = "fdr"))
)
adj <- 1 * (getmatrix(mod_train, "omega") != 0)
expect_equal(sum(adj[lower.tri(adj)]), 23)

mod_test_conf <- suppressWarnings(runmodel(ggm(test_data, omega = adj)))
expect_true(abs(mod_test_conf@fitmeasures$logl - (-26196.0963)) < 0.01)
expect_equal(mod_test_conf@fitmeasures$df, 22)
expect_true(abs(mod_test_conf@fitmeasures$chisq - 56.02942) < 0.01)

mod_test_sat <- suppressWarnings(runmodel(ggm(test_data)))
cmp2 <- compare(confirmatory = mod_test_conf, saturated = mod_test_sat)
expect_true(all(is.finite(cmp2$AIC)))
expect_true(all(is.finite(cmp2$BIC)))

## ---- Example 4: multi-group network comparison (gender) ----
bfig <- na.omit(bfi[, c(1:10, 26)])
bfig <- bfig[bfig$gender %in% c(1, 2), ]
mod_configural <- suppressWarnings(
  runmodel(ggm(bfig, vars = colnames(bfig)[1:10], groups = "gender"))
)
expect_true(abs(mod_configural@fitmeasures$logl - (-42265.3139)) < 0.01)
expect_equal(mod_configural@fitmeasures$npar, 130)
expect_equal(mod_configural@fitmeasures$df, 0)

mod_conf_pruned <- suppressWarnings(prune(mod_configural, alpha = 0.01))
expect_equal(mod_conf_pruned@fitmeasures$df, 48)

# Equality-constrained networks (applied to the consistent configural model):
mod_equal <- suppressWarnings(runmodel(groupequal(mod_configural, "omega")))
expect_equal(mod_equal@fitmeasures$df, 45)
cmp3 <- compare(configural = mod_configural, equal = mod_equal)
expect_true(all(is.finite(cmp3$AIC)))
expect_equal(cmp3$DF_diff[2], 45)
expect_true(is.finite(cmp3$p_value[2]))

omega_g1 <- getmatrix(mod_configural, "omega", group = 1)
omega_g2 <- getmatrix(mod_configural, "omega", group = 2)
expect_equal(dim(omega_g1), c(10L, 10L))
expect_equal(dim(omega_g2), c(10L, 10L))
expect_false(identical(omega_g1, omega_g2))

## ---- Example 7: GGM from a correlation matrix (psych::Thurstone) ----
mod_thu <- suppressWarnings(suppressMessages(
  runmodel(ggm(covs = psych::Thurstone, nobs = 213))
))
expect_true(abs(mod_thu@fitmeasures$logl - (-2166.56015)) < 0.01)
expect_equal(mod_thu@fitmeasures$df, 0)
expect_equal(mod_thu@fitmeasures$npar, 36)
mod_thu_pruned <- suppressWarnings(
  stepup(prune(mod_thu, alpha = 0.01, adjust = "fdr"))
)
expect_equal(mod_thu_pruned@fitmeasures$df, 23)
expect_equal(dim(getmatrix(mod_thu_pruned, "omega")), c(9L, 9L))

## ---- at_home: full 25-item versions of examples 1-2 ----
if (at_home()) {
  bfi25 <- na.omit(bfi[, 1:25])
  sat25 <- suppressWarnings(runmodel(ggm(bfi25)))
  expect_true(abs(sat25@fitmeasures$logl - (-97757.5045)) < 0.02)
  expect_equal(sat25@fitmeasures$npar, 350)
  pr25 <- suppressWarnings(prune(sat25, alpha = 0.01, adjust = "fdr"))
  expect_equal(pr25@fitmeasures$df, 195)
  expect_true(abs(pr25@fitmeasures$chisq - 597.1269) < 0.05)
  su25 <- suppressWarnings(stepup(pr25))
  expect_inherits(su25, "psychonetrics")
  expect_true(su25@fitmeasures$df <= pr25@fitmeasures$df)
  cmp25 <- compare(saturated = sat25, pruned = pr25, stepup = su25)
  expect_true(all(is.finite(cmp25$AIC)))
}
