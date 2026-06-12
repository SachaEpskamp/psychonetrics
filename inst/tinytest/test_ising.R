# Tests derived from psychonetrics.org examples 14 and 15 (Ising models on
# the bundled Jonas data, -1/1 coding), reduced to 8 (single group) and
# 6 (multi-group) nodes for speed. Also covers the multi-state Ising
# generalization (> 2 response options, added in 0.15.13).
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
data("Jonas", package = "psychonetrics")

## ---- Example 14: single-group Ising (8 items) ----
binary_items <- Jonas[, 1:8]
mod_sat <- suppressWarnings(runmodel(Ising(binary_items)))
expect_inherits(mod_sat, "psychonetrics")
expect_true(abs(mod_sat@fitmeasures$logl - (-773.88884)) < 0.01)
expect_equal(mod_sat@fitmeasures$npar, 36)  # 8 thresholds + 28 edges
expect_equal(mod_sat@fitmeasures$df, 0)

mod_pruned <- suppressWarnings(prune(mod_sat, alpha = 0.01, adjust = "fdr"))
expect_equal(mod_pruned@fitmeasures$df, 25)
expect_true(abs(mod_pruned@fitmeasures$chisq - 145.26824) < 0.01)

omega_ising <- getmatrix(mod_pruned, "omega")
expect_equal(dim(omega_ising), c(8L, 8L))
expect_true(all(diag(omega_ising) == 0))
tau <- getmatrix(mod_pruned, "tau")
expect_equal(length(as.vector(tau)), 8L)

cmp <- compare(saturated = mod_sat, pruned = mod_pruned)
expect_true(all(is.finite(cmp$AIC)))
expect_true(all(is.finite(cmp$BIC)))

## ---- Example 15 (reduced): multi-group Ising (6 items) ----
jit6 <- Jonas[, 1:6]
jit6$group <- Jonas$group
mod_mg <- suppressWarnings(
  runmodel(Ising(jit6, vars = colnames(jit6)[1:6], groups = "group"))
)
expect_true(abs(mod_mg@fitmeasures$logl - (-592.57710)) < 0.01)
expect_equal(mod_mg@fitmeasures$npar, 42)  # 2 x (6 tau + 15 omega)

mod_mg_eq <- suppressWarnings(runmodel(groupequal(mod_mg, "omega")))
expect_equal(mod_mg_eq@fitmeasures$df, 14)
cmp_mg <- compare(configural = mod_mg, equal_omega = mod_mg_eq)
expect_true(all(is.finite(cmp_mg$AIC)))
expect_equal(cmp_mg$DF_diff[2], 14)
expect_true(is.finite(cmp_mg$p_value[2]))

## ---- Multi-state Ising: 3 response options on simulated data ----
set.seed(42)
n <- 300
sim3 <- as.data.frame(matrix(
  sample(0:2, n * 4, replace = TRUE, prob = c(0.5, 0.3, 0.2)), n, 4
))
colnames(sim3) <- paste0("I", 1:4)
mod_ms <- suppressWarnings(runmodel(Ising(sim3)))
expect_inherits(mod_ms, "psychonetrics")
expect_true(abs(mod_ms@fitmeasures$logl - (-1217.91247)) < 0.01)
expect_equal(mod_ms@fitmeasures$npar, 10)  # 4 thresholds + 6 edges
expect_equal(mod_ms@fitmeasures$df, 0)

## ---- at_home: example 14/15 fuller (10 items + stepup; equal thresholds) ----
if (at_home()) {
  mod_sat10 <- suppressWarnings(runmodel(Ising(Jonas[, 1:10])))
  expect_true(is.finite(mod_sat10@fitmeasures$logl))
  mod_pr10 <- suppressWarnings(prune(mod_sat10, alpha = 0.01, adjust = "fdr"))
  mod_su10 <- suppressWarnings(stepup(mod_pr10))
  expect_true(mod_su10@fitmeasures$df <= mod_pr10@fitmeasures$df)
  cmp10 <- compare(saturated = mod_sat10, pruned = mod_pr10, stepup = mod_su10)
  expect_true(all(is.finite(cmp10$AIC)))

  # Further equality constraints on the multi-group model:
  mod_mg_eq2 <- suppressWarnings(runmodel(groupequal(mod_mg_eq, "tau")))
  expect_true(mod_mg_eq2@fitmeasures$df > mod_mg_eq@fitmeasures$df)
  expect_true(is.finite(mod_mg_eq2@fitmeasures$logl))
}
