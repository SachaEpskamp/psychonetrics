# Tests for the meta-analytic varcov family (meta_varcov / MAGNA):
# pooled correlation recovery, the bivariate (2x2) case, estimator = "ML",
# and Vmethod = "OS_individual". Data are simulated with fixed seeds.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))

## ---- Pooled correlation recovery, 8 studies x 4 variables ----
set.seed(123)
nv <- 4; nstudy <- 8
trueSigma <- matrix(0.3, nv, nv); diag(trueSigma) <- 1
corList <- vector("list", nstudy); ns <- numeric(nstudy)
for (s in 1:nstudy) {
  ns[s] <- sample(150:400, 1)
  X <- MASS::mvrnorm(ns[s], rep(0, nv), trueSigma)
  corList[[s]] <- cor(X)
}

mod_meta <- suppressWarnings(runmodel(meta_varcov(cors = corList, nobs = ns)))
expect_inherits(mod_meta, "psychonetrics")
expect_true(abs(mod_meta@fitmeasures$logl - 76.03953) < 0.01)
expect_equal(mod_meta@fitmeasures$npar, 27)
expect_equal(mod_meta@fitmeasures$df, 0)

rho <- getmatrix(mod_meta, "rho_y")
expect_equal(dim(rho), c(4L, 4L))
# Pooled correlations recover the generating value (0.3):
expect_true(max(abs(rho[lower.tri(rho)] - 0.3)) < 0.06)

## ---- Bivariate (2x2) meta_varcov (regression: must run, recover r) ----
set.seed(99)
corList2 <- vector("list", 6); ns2 <- numeric(6)
for (s in 1:6) {
  ns2[s] <- sample(100:300, 1)
  X <- MASS::mvrnorm(ns2[s], c(0, 0), matrix(c(1, 0.4, 0.4, 1), 2, 2))
  corList2[[s]] <- cor(X)
}
mod_meta2 <- suppressWarnings(runmodel(meta_varcov(cors = corList2, nobs = ns2)))
expect_inherits(mod_meta2, "psychonetrics")
expect_true(abs(mod_meta2@fitmeasures$logl - 8.55549) < 0.01)
expect_true(abs(getmatrix(mod_meta2, "rho_y")[2, 1] - 0.40252) < 0.01)

## ---- estimator = "ML" (regression: must run; equals FIML w/o missings) ----
mod_meta_ml <- suppressWarnings(
  runmodel(meta_varcov(cors = corList2, nobs = ns2, estimator = "ML"))
)
expect_inherits(mod_meta_ml, "psychonetrics")
expect_true(abs(mod_meta_ml@fitmeasures$logl - 8.55549) < 0.01)

## ---- Vmethod = "OS_individual" (regression: must run) ----
mod_meta_os <- suppressWarnings(
  runmodel(meta_varcov(cors = corList2, nobs = ns2, Vmethod = "OS_individual"))
)
expect_inherits(mod_meta_os, "psychonetrics")
expect_true(abs(mod_meta_os@fitmeasures$logl - 8.55057) < 0.01)
