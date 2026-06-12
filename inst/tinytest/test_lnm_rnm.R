# Tests derived from psychonetrics.org examples 9, 10 and 17 (latent and
# residual network models on the bundled StarWars data).
#
# Note: lnm() / rnm() already hardcode latent = "ggm" / residual = "ggm";
# the psychonetrics.org examples pass those arguments again, which errors
# ("matched by multiple actual arguments") - they are omitted here.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
data("StarWars", package = "psychonetrics")

items <- StarWars[, 1:7]
lambda <- matrix(0, 7, 3)
lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1; lambda[7, 3] <- 1
colnames(lambda) <- c("Original", "Prequel", "Sequel")

## ---- Example 9: latent network model (LNM) ----
mod_lnm <- suppressWarnings(
  runmodel(lnm(items, lambda = lambda, identification = "variance"))
)
expect_inherits(mod_lnm, "psychonetrics")
expect_true(abs(mod_lnm@fitmeasures$logl - (-2939.63901)) < 0.01)
expect_equal(mod_lnm@fitmeasures$df, 12)
expect_equal(mod_lnm@fitmeasures$npar, 23)

mod_lnm_pruned <- suppressWarnings(
  stepup(prune(mod_lnm, alpha = 0.01, adjust = "fdr", matrices = "omega_zeta"),
         matrices = "omega_zeta")
)
expect_equal(mod_lnm_pruned@fitmeasures$df, 13)

cmp <- compare(saturated = mod_lnm, pruned = mod_lnm_pruned)
expect_true(all(is.finite(cmp$AIC)))
expect_true(all(is.finite(cmp$BIC)))

omega_latent <- getmatrix(mod_lnm_pruned, "omega_zeta")
expect_equal(dim(omega_latent), c(3L, 3L))
expect_true(all(diag(omega_latent) == 0))

## ---- Example 10: residual network model (RNM) ----
# omega_epsilon defaults to "zero" (an empty residual network):
mod_rnm_empty <- suppressWarnings(
  runmodel(rnm(items, lambda = lambda, omega_epsilon = "zero",
               identification = "variance"))
)
expect_true(abs(mod_rnm_empty@fitmeasures$logl - (-2939.63901)) < 0.01)
expect_equal(mod_rnm_empty@fitmeasures$df, 12)

mod_rnm_stepup <- suppressWarnings(
  stepup(mod_rnm_empty, matrices = "omega_epsilon", alpha = 0.01)
)
expect_equal(mod_rnm_stepup@fitmeasures$df, 11)
omega_residual <- getmatrix(mod_rnm_stepup, "omega_epsilon")
expect_equal(dim(omega_residual), c(7L, 7L))
expect_true(all(diag(omega_residual) == 0))

# Comparison with a standard CFA:
mod_cfa7 <- suppressWarnings(
  runmodel(lvm(items, lambda = lambda, identification = "variance"))
)
cmp2 <- compare(standard_CFA = mod_cfa7, residual_network = mod_rnm_stepup)
expect_true(all(is.finite(cmp2$AIC)))

## ---- Example 17 (reduced): CFA on 10 items + freepar ----
obsvars <- paste0("Q", 1:10)
latents <- c("Prequels", "Original", "Sequels")
Lambda <- matrix(0, 10, 3)
Lambda[1:4, 1] <- 1; Lambda[c(1, 5:7), 2] <- 1; Lambda[c(1, 8:10), 3] <- 1

mod1 <- suppressWarnings(
  runmodel(lvm(StarWars, lambda = Lambda, vars = obsvars,
               identification = "variance", latents = latents))
)
expect_true(abs(mod1@fitmeasures$logl - (-4126.35127)) < 0.01)
expect_equal(mod1@fitmeasures$df, 30)
expect_true(abs(mod1@fitmeasures$chisq - 34.57912) < 0.01)

# Free the Q10-Q4 residual covariance (one extra parameter, lower chisq):
mod2 <- suppressWarnings(
  runmodel(freepar(mod1, "sigma_epsilon", "Q10", "Q4"))
)
expect_equal(mod2@fitmeasures$df, 29)
expect_true(abs(mod2@fitmeasures$chisq - 25.07150) < 0.01)
expect_true(mod2@fitmeasures$chisq < mod1@fitmeasures$chisq)

cmp3 <- compare(original = mod1, adjusted = mod2)
expect_equal(abs(cmp3$DF[1] - cmp3$DF[2]), 1)
expect_true(all(is.finite(cmp3$AIC)))
expect_true(is.finite(cmp3$p_value[2]))

## ---- at_home: example 17 continued (LNM on 10 items + prune) ----
if (at_home()) {
  lnm_mod <- suppressWarnings(
    prune(runmodel(lnm(StarWars, lambda = Lambda, vars = obsvars,
                       latents = latents, identification = "variance")),
          alpha = 0.05)
  )
  expect_inherits(lnm_mod, "psychonetrics")
  expect_true(is.finite(lnm_mod@fitmeasures$logl))
  cmp4 <- compare(cfa = mod1, lnm = lnm_mod)
  expect_true(all(is.finite(cmp4$AIC)))

  # Exploratory GGM from summary statistics (cov matrix + nobs):
  S <- (nrow(StarWars) - 1) / nrow(StarWars) * cov(StarWars[, 1:10])
  ggm_mod <- suppressWarnings(suppressMessages(
    prune(stepup(ggm(covs = S, nobs = nrow(StarWars), omega = "zero",
                     covtype = "ML")))
  ))
  expect_inherits(ggm_mod, "psychonetrics")
  expect_true(is.finite(ggm_mod@fitmeasures$logl))
}
