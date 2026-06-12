# Tests derived from psychonetrics.org examples 6 and 16 (FIML estimation),
# reduced to 10 bfi items. Also covers the per-row fullFIML = TRUE path,
# which must match the missingness-pattern FIML likelihood exactly.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
if (!requireNamespace("psych", quietly = TRUE)) exit_file("psych not available")

bfi <- psych::bfi

## ---- Example 6: GGM with missing data (FIML via missing = "auto") ----
bfi10m <- bfi[, 1:10]
expect_true(anyNA(bfi10m))

mod_fiml <- suppressWarnings(suppressMessages(runmodel(ggm(bfi10m, missing = "auto"))))
expect_inherits(mod_fiml, "psychonetrics")
expect_equal(mod_fiml@estimator, "FIML")
expect_true(abs(mod_fiml@fitmeasures$logl - (-44796.3392)) < 0.01)
expect_equal(mod_fiml@fitmeasures$npar, 65)

mod_fiml_pruned <- suppressWarnings(prune(mod_fiml, alpha = 0.01, adjust = "fdr"))
expect_equal(mod_fiml_pruned@fitmeasures$df, 21)
omega_fiml <- getmatrix(mod_fiml_pruned, "omega")
expect_equal(dim(omega_fiml), c(10L, 10L))

## ---- fullFIML = TRUE (per-row likelihood) equals pattern-based FIML ----
tiny <- bfi[1:300, 1:4]
m_pattern <- suppressWarnings(runmodel(ggm(tiny, estimator = "FIML")))
m_perrow  <- suppressWarnings(runmodel(ggm(tiny, estimator = "FIML", fullFIML = TRUE)))
expect_true(abs(m_pattern@fitmeasures$logl - (-1908.21589)) < 0.01)
expect_true(abs(m_pattern@fitmeasures$logl - m_perrow@fitmeasures$logl) < 1e-6)
expect_equal(m_pattern@fitmeasures$npar, m_perrow@fitmeasures$npar)

## ---- at_home: example 16 (NA2020 sleep & wellbeing, FIML) ----
if (at_home()) {
  data("NA2020", package = "psychonetrics")
  mod <- suppressWarnings(runmodel(ggm(NA2020, estimator = "FIML")))
  expect_true(abs(mod@fitmeasures$logl - (-7392.44263)) < 0.01)

  network <- getmatrix(mod, "omega", threshold = TRUE, alpha = 0.05)
  expect_equal(dim(network), c(8L, 8L))

  mod_prune <- suppressWarnings(prune(mod, alpha = 0.05))
  expect_equal(mod_prune@fitmeasures$df, 16)
  mod_stepup <- suppressWarnings(stepup(mod_prune))
  expect_equal(mod_stepup@fitmeasures$df, 15)

  mod_empty_stepup <- suppressWarnings(
    stepup(runmodel(ggm(NA2020, estimator = "FIML", omega = "zero")))
  )
  expect_inherits(mod_empty_stepup, "psychonetrics")
  expect_true(is.finite(mod_empty_stepup@fitmeasures$logl))

  # (suppressWarnings: the from_empty stepup model is not nested in the others,
  # so compare() warns about a negative chi-square difference)
  cmp <- suppressWarnings(
    compare(saturated = mod, pruned = mod_prune, stepup = mod_stepup,
            from_empty = mod_empty_stepup)
  )
  expect_true(all(is.finite(cmp$AIC)))
  expect_true(all(is.finite(cmp$BIC)))
}
