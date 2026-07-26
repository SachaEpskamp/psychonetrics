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

## ---- (e) The verdict table (0.17.2 guarded decision tree) ----
# Default criterion is the chi-square difference test (since 0.17.5):
expect_equal(eval(formals(ri_clpm_search)$criterion)[1], "Chisq")
expect_equal(suppressWarnings(ri_clpm_search(fit_cor, verbose = FALSE))$criterion, "Chisq")
expect_true(is.data.frame(res$verdicts))
expect_equal(res$verdicts$step,
             c("temporal", "contemporaneous", "intercepts", "innovation",
               "wave1_endogenous"))
expect_true(all(c("step", "tested", "decision", "Chisq_diff", "DF_diff",
                  "p_value", "dAIC", "dBIC", "note") %in% names(res$verdicts)))
expect_true(all(res$verdicts$decision %in%
                  c("retained", "released", "rejected (gate)", "not tested",
                    "not applicable")))
expect_true(is.character(res$selected_name))
expect_equal(res$on_rejection, "release")

# Every step after a retained temporal gate must have been tested (release
# rule: only the gate can stop the search). The wave-1 step may be "not
# applicable" when the covariance process is not stationary:
if (res$verdicts$decision[1] == "retained"){
  expect_true(all(res$verdicts$tested[2:4]))
} else {
  expect_equal(res$verdicts$decision[1], "rejected (gate)")
  expect_true(all(res$verdicts$decision[2:5] == "not tested"))
}

# The wave-1 step is gated on covariance stationarity:
covstat <- all(res$verdicts$decision[c(1, 2, 4)] == "retained")
if (!covstat){
  expect_true(res$verdicts$decision[5] %in% c("not applicable", "not tested"))
} else if (isTRUE(res$verdicts$tested[5])){
  # The panel model is the CONSTRAINED special case: m(m+1)/2 more df:
  expect_equal(res$verdicts$DF_diff[5], nVar * (nVar + 1) / 2)
  expect_true(res$verdicts$Chisq_diff[5] >= -1e-4)
}

## ---- (f) on_rejection = "stop" reproduces the old chain ----
res_stop <- suppressWarnings(ri_clpm_search(fit_cor, criterion = "Chisq",
                                            verbose = FALSE,
                                            on_rejection = "stop"))
expect_equal(res_stop$on_rejection, "stop")
# In "stop" mode no step after a rejection is ever tested:
dec <- res_stop$verdicts$decision
firstreject <- which(dec %in% c("released", "rejected (gate)"))
if (length(firstreject) > 0){
  expect_true(all(dec[(firstreject[1] + 1):length(dec)] == "not tested"))
}

## ---- (g) RI-CLPN wrappers ----
expect_equal(eval(formals(ri_clpn)$type)[1], "ggm")
mod_clpn <- suppressWarnings(ri_clpn(dat, design))
mod_ggm  <- suppressWarnings(ri_clpm(dat, design, type = "ggm"))
expect_equal(mod_clpn@extramatrices$ri_clpm_type, "ggm")
expect_identical(mod_clpn@parameters, mod_ggm@parameters)
# The type argument can still be overridden:
expect_equal(suppressWarnings(ri_clpn(dat, design, type = "cov"))@extramatrices$ri_clpm_type,
             "cov")
# ri_clpn_stationary()/ri_clpn_search() are pure aliases:
expect_identical(suppressWarnings(ri_clpn_stationary(fit_cor, "temporal"))@parameters,
                 suppressWarnings(ri_clpm_stationary(fit_cor, "temporal"))@parameters)
res_alias <- suppressWarnings(ri_clpn_search(fit_cor, criterion = "Chisq", verbose = FALSE))
expect_identical(res_alias$verdicts, res$verdicts)
expect_identical(res_alias$selected_name, res$selected_name)

## ---- (h) Block extraction: getmatrix(..., time = ) (0.17.3) ----
nObs <- nVar * nWave
lat  <- function(t) if (is.na(t)) nObs + seq_len(nVar) else (t - 1) * nVar + seq_len(nVar)

full_beta <- getmatrix(fit_cor, "beta")
full_rho  <- getmatrix(fit_cor, "rho_zeta")

# Plain names without 'time' are unchanged (backward compatibility):
expect_equal(dim(full_beta), c(nObs + nVar, nObs + nVar))

# Temporal blocks:
for (tt in seq_len(nWave - 1)){
  expect_equal(unname(getmatrix(fit_cor, "beta", time = tt)),
               unname(full_beta[lat(tt + 1), lat(tt)]))
}
# Contemporaneous and between-person blocks:
expect_equal(unname(getmatrix(fit_cor, "rho_zeta_within", time = 2)),
             unname(full_rho[lat(2), lat(2)]))
expect_equal(unname(getmatrix(fit_cor, "rho_zeta_between")),
             unname(full_rho[lat(NA), lat(NA)]))
# Shorthand and wave labels:
expect_equal(getmatrix(fit_cor, "rho_zeta", time = 2),
             getmatrix(fit_cor, "rho_zeta_within", time = 2))
expect_equal(getmatrix(fit_cor, "rho_zeta_within", time = "T3"),
             getmatrix(fit_cor, "rho_zeta_within", time = 3))
# Blocks carry the design matrix variable names:
expect_equal(rownames(getmatrix(fit_cor, "rho_zeta_between")), paste0("V", 1:nVar))
expect_equal(colnames(getmatrix(fit_cor, "beta", time = 1)), paste0("V", 1:nVar))
# Mean structure:
expect_equal(unname(getmatrix(fit_cor, "nu", time = 2)),
             unname(getmatrix(fit_cor, "nu")[lat(2), , drop = FALSE]))
expect_equal(unname(getmatrix(fit_cor, "nu_between")),
             unname(getmatrix(fit_cor, "nu_eta")[lat(NA), , drop = FALSE]))

## ---- (i) PDC blocks ----
manualPDC <- function(tt){
  bb <- full_beta[lat(tt + 1), lat(tt)]
  kk <- getmatrix(fit_cor, "kappa_zeta")[lat(tt + 1), lat(tt + 1)]
  unname(t(bb / sqrt(diag(solve(kk)) %o% diag(kk) + bb^2)))
}
for (tt in seq_len(nWave - 1)){
  expect_equal(unname(getmatrix(fit_cor, "PDC", time = tt)), manualPDC(tt))
}
expect_equal(dimnames(getmatrix(fit_cor, "PDC", time = 1)),
             list(paste0("V", 1:nVar), paste0("V", 1:nVar)))

## ---- (j) Thresholding is exactly the submatrix of the thresholded matrix ----
thr_full  <- suppressWarnings(getmatrix(fit_cor, "rho_zeta", threshold = TRUE, alpha = 0.05))
thr_block <- suppressWarnings(getmatrix(fit_cor, "rho_zeta_within", time = 2,
                                        threshold = TRUE, alpha = 0.05))
expect_equal(unname(thr_block), unname(thr_full[lat(2), lat(2)]))
# Numeric thresholds:
raw2 <- getmatrix(fit_cor, "rho_zeta_within", time = 2)
ref2 <- raw2; ref2[abs(ref2) < 0.05] <- 0
expect_equal(getmatrix(fit_cor, "rho_zeta_within", time = 2, threshold = 0.05), ref2)

## ---- (k) diag = FALSE applies to the block ----
b1  <- getmatrix(fit_cor, "beta", time = 1)
b1d <- getmatrix(fit_cor, "beta", time = 1, diag = FALSE)
expect_true(all(diag(b1d) == 0))
expect_true(any(diag(b1) != 0))
expect_equal(b1d[lower.tri(b1d)], b1[lower.tri(b1)])

## ---- (l) Collapsing over waves after a stationarity constraint ----
stat_t <- suppressWarnings(runmodel(ri_clpm_stationary(fit_cor, "temporal"), verbose = FALSE))
expect_equal(getmatrix(stat_t, "beta", time = 1), getmatrix(stat_t, "beta", time = 2))
# Unconstrained: one block per wave; the between block is a single matrix:
expect_true(is.list(getmatrix(fit_cor, "PDC")))
expect_equal(length(getmatrix(fit_cor, "PDC")), nWave - 1L)
expect_true(is.matrix(getmatrix(fit_cor, "rho_zeta_between")))

## ---- (m) Guards ----
expect_error(getmatrix(fit_cor, "rho_zeta_between", time = 1))
expect_error(getmatrix(fit_cor, "beta", time = nWave))
expect_error(getmatrix(fit_cor, "rho_zeta_within", time = nWave + 1))
expect_error(getmatrix(fit_cor, "rho_zeta_within", time = "T99"))
