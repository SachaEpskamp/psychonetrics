# Regression tests for the RI-CLPM helpers ri_clpm(), ri_clpm_stationary()
# and ri_clpm_search(). Covers (a) the default type = "cor" parameterization,
# (b) that the unconstrained "cor" and "cov" models are equivalent in fit,
# (c) that the innovation-variance and contemporaneous-correlation
# stationarity constraints are cleanly separated and never touch the
# exogenous wave-1 block nor the between-person random-intercept block, and
# (d) that ri_clpm_search() visits the constraints in the intended order
# (temporal -> contemporaneous -> intercepts -> innovation -> panel(G)VAR).
#
# Data are simulated with a fixed seed (3 variables, 4 waves).

suppressMessages(library(psychonetrics))

## ---- Simulate a small random-intercept panel data set ----
set.seed(1)
nPerson <- 250; nVar <- 3; nWave <- 4
RI <- matrix(rnorm(nPerson * nVar), nPerson, nVar)
beta <- diag(0.3, nVar); beta[2, 1] <- 0.2; beta[3, 2] <- 0.2
within <- array(0, c(nPerson, nVar, nWave))
within[, , 1] <- matrix(rnorm(nPerson * nVar), nPerson, nVar)
for (w in 2:nWave){
  within[, , w] <- within[, , w - 1] %*% t(beta) +
    matrix(rnorm(nPerson * nVar, sd = 0.8), nPerson, nVar)
}
dat <- do.call(cbind, lapply(1:nWave, function(w) within[, , w] + RI))
colnames(dat) <- as.vector(outer(paste0("V", 1:nVar), 1:nWave, paste, sep = "_"))
dat <- as.data.frame(dat)
design <- matrix(colnames(dat), nVar, nWave)
rownames(design) <- paste0("V", 1:nVar)

n_obs <- nVar * nWave

## ---- (a) Default type is "cor" ----
expect_equal(eval(formals(ri_clpm)$type)[1], "cor")
mod_cor <- suppressWarnings(ri_clpm(dat, design))
expect_equal(mod_cor@extramatrices$ri_clpm_type, "cor")

## ---- (b) Unconstrained "cor" and "cov" have identical fit ----
mod_cov <- suppressWarnings(ri_clpm(dat, design, type = "cov"))
fit_cor <- suppressWarnings(runmodel(mod_cor, verbose = FALSE))
fit_cov <- suppressWarnings(runmodel(mod_cov, verbose = FALSE))
expect_true(abs(fit_cor@fitmeasures$chisq - fit_cov@fitmeasures$chisq) < 1e-3)
expect_equal(fit_cor@fitmeasures$df, fit_cov@fitmeasures$df)
expect_equal(max(fit_cor@parameters$par), max(fit_cov@parameters$par))

## ---- (c) "cor" splits into rho_zeta (off-diagonal) + SD_zeta (diagonal) ----
pp  <- fit_cor@parameters
zet <- pp[grepl("_zeta", pp$matrix) & !pp$fixed, ]
expect_equal(sort(unique(zet$matrix)), c("rho_zeta", "SD_zeta"))
expect_equal(unique(zet$matrix[zet$row != zet$col]), "rho_zeta")
expect_equal(unique(zet$matrix[zet$row == zet$col]), "SD_zeta")

## Helper: block index of a latent (1..nWave = wave blocks, nWave+1 = RI block)
blockOf <- function(idx) ifelse(idx > n_obs, nWave + 1L, ceiling(idx / nVar))

## ---- (d) Innovation-variance constraint ----
inno <- suppressWarnings(ri_clpm_stationary(fit_cor, "innovation"))
si <- inno@parameters
si <- si[si$matrix == "SD_zeta" & !si$fixed & si$row == si$col, ]
si$block <- blockOf(si$row)
inno_pars <- unique(si$par[si$block >= 2 & si$block <= nWave])   # innovation blocks
w1_pars   <- unique(si$par[si$block == 1])                        # exogenous wave 1
ri_pars   <- unique(si$par[si$block == nWave + 1])                # between-person RI
# Innovation variances are equated across waves (one par per variable) ...
si$within <- si$row - ifelse(si$block > nWave, nWave * nVar, (si$block - 1) * nVar)
for (v in seq_len(nVar)){
  pv <- si$par[si$within == v & si$block >= 2 & si$block <= nWave]
  expect_equal(length(unique(pv)), 1L)
}
# ... but never equated with the wave-1 exogenous or the RI variances:
expect_equal(length(intersect(inno_pars, w1_pars)), 0L)
expect_equal(length(intersect(inno_pars, ri_pars)), 0L)
# Degrees of freedom added = nVar * (nWave - 2):
inno_fit <- suppressWarnings(runmodel(inno, verbose = FALSE))
expect_equal(inno_fit@fitmeasures$df - fit_cor@fitmeasures$df, nVar * (nWave - 2))

## ---- (d) Contemporaneous-correlation constraint ----
con <- suppressWarnings(ri_clpm_stationary(fit_cor, "contemporaneous"))
con_fit <- suppressWarnings(runmodel(con, verbose = FALSE))
addexp <- (nVar * (nVar - 1) / 2) * (nWave - 2)
expect_equal(con_fit@fitmeasures$df - fit_cor@fitmeasures$df, addexp)

## ---- (d) ri_clpm_search visits the constraints in the intended order ----
res <- suppressWarnings(ri_clpm_search(fit_cor, criterion = "Chisq", verbose = FALSE))
expect_inherits(res, "ri_clpm_search")
constrained <- res$path$model[res$path$model != "base"]
expected <- c("temporal", "contemporaneous", "intercepts", "innovation", "panelVAR")
# Whatever steps were attempted must be a prefix of the intended order:
expect_true(length(constrained) >= 1L)
expect_equal(constrained, expected[seq_along(constrained)])
