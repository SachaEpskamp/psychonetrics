# Tests derived from psychonetrics.org examples 8 (CFA) and 11 (SEM).
#
# Note on example 11: the psychonetrics.org code sets diag(sigma_epsilon) <- 0,
# which FIXES all residual variances to zero (and frees the full y1-y4
# residual covariance block), yielding a non-identified model that does not
# converge. These tests use the textbook PoliticalDemocracy specification
# (free residual variances + the six classic residual covariances), which
# psychonetrics fits identically to lavaan.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")

## ---- Example 8: CFA on HolzingerSwineford1939 ----
hs_data <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
lambda <- matrix(0, 9, 3)
lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1; lambda[7:9, 3] <- 1
colnames(lambda) <- c("visual", "textual", "speed")

mod_cfa <- suppressWarnings(
  runmodel(lvm(hs_data, lambda = lambda, identification = "variance"))
)
expect_inherits(mod_cfa, "psychonetrics")
expect_true(abs(mod_cfa@fitmeasures$logl - (-3737.74493)) < 0.01)
expect_true(abs(mod_cfa@fitmeasures$chisq - 85.30553) < 0.01)
expect_equal(mod_cfa@fitmeasures$df, 24)
expect_equal(mod_cfa@fitmeasures$npar, 30)
expect_true(abs(mod_cfa@fitmeasures$rmsea - 0.092121) < 1e-3)
expect_true(abs(mod_cfa@fitmeasures$cfi - 0.930560) < 1e-3)

# Cross-check the log-likelihood against lavaan (computed in-test):
lav_fit <- lavaan::cfa(
  "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6\nspeed =~ x7 + x8 + x9",
  data = hs_data, std.lv = TRUE, meanstructure = TRUE
)
lav_ll <- unname(lavaan::fitMeasures(lav_fit)["logl"])
expect_true(abs(mod_cfa@fitmeasures$logl - lav_ll) < 1e-3)

# Latent correlation matrix: variance identification fixes the diagonal to 1:
sigma_zeta <- getmatrix(mod_cfa, "sigma_zeta")
expect_equal(dim(sigma_zeta), c(3L, 3L))
expect_true(all(abs(diag(sigma_zeta) - 1) < 1e-8))
expect_true(isSymmetric(unname(sigma_zeta)))

# Residual variances are positive:
residual_var <- diag(getmatrix(mod_cfa, "sigma_epsilon"))
expect_equal(length(residual_var), 9L)
expect_true(all(residual_var > 0))

# MIs returns a data frame invisibly:
invisible(capture.output(mi <- MIs(mod_cfa, verbose = FALSE)))
expect_true(is.data.frame(mi))
expect_true(nrow(mi) > 0)

## ---- Example 11: SEM with structural paths (PoliticalDemocracy) ----
vars <- c(paste0("y", 1:8), paste0("x", 1:3))
data_sem <- lavaan::PoliticalDemocracy[, vars]
lambda2 <- matrix(0, 11, 3)
lambda2[9:11, 1] <- 1   # ind60 =~ x1 + x2 + x3
lambda2[1:4, 2]  <- 1   # dem60 =~ y1 + y2 + y3 + y4
lambda2[5:8, 3]  <- 1   # dem65 =~ y5 + y6 + y7 + y8
colnames(lambda2) <- c("ind60", "dem60", "dem65")
rownames(lambda2) <- vars
beta <- matrix(0, 3, 3)
beta[2, 1] <- 1          # dem60 ~ ind60
beta[3, 1] <- 1          # dem65 ~ ind60
beta[3, 2] <- 1          # dem65 ~ dem60
rownames(beta) <- colnames(beta) <- c("ind60", "dem60", "dem65")
# Textbook correlated residuals: y1-y5, y2-y4, y2-y6, y3-y7, y4-y8, y6-y8.
# Residual variances (the diagonal) must be free!
sigma_epsilon <- matrix(0, 11, 11)
res_pairs <- rbind(c(1, 5), c(2, 4), c(2, 6), c(3, 7), c(4, 8), c(6, 8))
for (i in seq_len(nrow(res_pairs))) {
  sigma_epsilon[res_pairs[i, 1], res_pairs[i, 2]] <- 1
  sigma_epsilon[res_pairs[i, 2], res_pairs[i, 1]] <- 1
}
diag(sigma_epsilon) <- 1

mod_sem <- suppressWarnings(
  runmodel(lvm(data_sem, lambda = lambda2, beta = beta,
               sigma_epsilon = sigma_epsilon,
               identification = "variance", vars = vars))
)
expect_inherits(mod_sem, "psychonetrics")
# This model previously diverged; it must now converge to the lavaan solution:
expect_true(is.finite(mod_sem@fitmeasures$logl))
expect_true(abs(mod_sem@fitmeasures$logl - (-1547.79096)) < 0.01)
expect_true(abs(mod_sem@fitmeasures$chisq - 38.12524) < 0.01)
expect_equal(mod_sem@fitmeasures$df, 35)
expect_equal(mod_sem@fitmeasures$npar, 42)

# Cross-check against lavaan:
lav_sem <- lavaan::sem(
  "ind60 =~ x1 + x2 + x3
   dem60 =~ y1 + y2 + y3 + y4
   dem65 =~ y5 + y6 + y7 + y8
   dem60 ~ ind60
   dem65 ~ ind60 + dem60
   y1 ~~ y5
   y2 ~~ y4 + y6
   y3 ~~ y7
   y4 ~~ y8
   y6 ~~ y8",
  data = data_sem, std.lv = TRUE, meanstructure = TRUE
)
lav_sem_ll <- unname(lavaan::fitMeasures(lav_sem)["logl"])
expect_true(abs(mod_sem@fitmeasures$logl - lav_sem_ll) < 1e-3)

sigma_zeta2 <- getmatrix(mod_sem, "sigma_zeta")
expect_equal(dim(sigma_zeta2), c(3L, 3L))
invisible(capture.output(pars_sem <- parameters(mod_sem)))
expect_true(is.data.frame(pars_sem))
