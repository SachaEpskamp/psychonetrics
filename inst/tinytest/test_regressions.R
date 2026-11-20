# Targeted regression tests for bugs fixed during the 0.14.x-0.15.x audit
# series. Each block is self-contained and deterministic.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
if (!requireNamespace("psych", quietly = TRUE)) exit_file("psych not available")

bfi <- psych::bfi

## ---- covtype = "choose" detects the UB denominator for cov(X) input ----
# (Integer-valued data allow the denominator detection to succeed.)
bfi6 <- na.omit(bfi[, 1:6])
S <- cov(bfi6)
m_choose <- suppressWarnings(suppressMessages(
  runmodel(ggm(covs = S, nobs = nrow(bfi6), covtype = "choose"))
))
m_ub <- suppressWarnings(suppressMessages(
  runmodel(ggm(covs = S, nobs = nrow(bfi6), covtype = "UB"))
))
m_ml <- suppressWarnings(suppressMessages(
  runmodel(ggm(covs = S, nobs = nrow(bfi6), covtype = "ML"))
))
expect_true(abs(m_choose@fitmeasures$logl - m_ub@fitmeasures$logl) < 1e-8)
expect_true(abs(m_choose@fitmeasures$logl - m_ml@fitmeasures$logl) > 0.1)

## ---- getmatrix(..., diag = FALSE) zeroes the diagonal ----
sg <- getmatrix(m_choose, "sigma", diag = FALSE)
expect_true(all(diag(sg) == 0))
sg2 <- getmatrix(m_choose, "sigma")  # default keeps the diagonal
expect_true(all(diag(sg2) > 0))

## ---- ebic.75 follows the EBIC formula with gamma = 0.75 ----
fm <- m_choose@fitmeasures
n_total <- sum(m_choose@sample@groups$nobs)
manual_ebic75 <- -2 * fm$logl + fm$npar * log(n_total) +
  4 * fm$npar * 0.75 * log(fm$nvar)
expect_true(abs(fm$ebic.75 - manual_ebic75) < 1e-6)

## ---- compare() reports NA p-value for an equal-df comparison ----
cmp_eq <- compare(a = m_choose, b = m_choose)
expect_equal(cmp_eq$DF_diff[2], 0)
expect_true(is.na(cmp_eq$p_value[2]))

## ---- spectralshift: min eigenvalue 0.001, off-diagonals untouched ----
A <- matrix(c(1, 0.9, 0.9,  0.9, 1, 0.9,  0.9, 0.9, -0.5), 3, 3)
expect_true(min(eigen(A)$values) < 0)  # indefinite input
As <- psychonetrics:::spectralshift(A)
expect_true(abs(min(eigen(As)$values) - 0.001) < 1e-10)
expect_equal(As[lower.tri(As)], A[lower.tri(A)])

## ---- prune(adjust = "bonferroni") adjusts over UNIQUE parameters ----
# Two-group model with omega equal across groups: each edge is one unique
# parameter; the bonferroni correction must count it once, not per group.
bfig <- na.omit(bfi[, c(1:10, 26)])
bfig <- bfig[bfig$gender %in% c(1, 2), ]
m_eq <- suppressWarnings(runmodel(groupequal(
  ggm(bfig, vars = colnames(bfig)[1:10], groups = "gender"), "omega"
)))
pt <- m_eq@parameters
expect_true(is.data.frame(pt))  # parameter table stays a data.frame
om <- pt[pt$matrix == "omega" & !pt$fixed, ]
p_unique <- om$p[!duplicated(om$par)]
expected_kept <- sum(p.adjust(p_unique, "bonferroni") < 0.01)
m_pruned <- suppressWarnings(prune(m_eq, alpha = 0.01, adjust = "bonferroni"))
pt2 <- m_pruned@parameters
om2 <- pt2[pt2$matrix == "omega" & !pt2$fixed, ]
expect_equal(sum(!duplicated(om2$par)), expected_kept)

## ---- parameters tables stay data.frame after groupequal() ----
expect_true(is.data.frame(m_eq@parameters))
invisible(capture.output(pars_eq <- parameters(m_eq)))
expect_true(is.data.frame(pars_eq))

## ---- equal = "omega_zeta" in a 2-group lnm shares parameter indices ----
data("StarWars", package = "psychonetrics")
sw <- StarWars[, 1:7]
sw$agegroup <- ifelse(StarWars$Q12 < 30, "young", "older")
lambda7 <- matrix(0, 7, 3)
lambda7[1:3, 1] <- 1; lambda7[4:6, 2] <- 1; lambda7[7, 3] <- 1
m_lnm_eq <- suppressWarnings(
  lnm(sw, vars = paste0("Q", 1:7), lambda = lambda7, groups = "agegroup",
      identification = "variance", equal = "omega_zeta")
)
oz <- m_lnm_eq@parameters
oz <- oz[oz$matrix == "omega_zeta" & oz$par > 0, ]
expect_true(nrow(oz) > 0)
shared <- vapply(
  split(oz$par, paste(oz$row, oz$col)),
  function(z) length(unique(z)) == 1L, logical(1)
)
expect_true(all(shared))

## ---- factorscores with beta != 0 matches the manual regression formula ----
set.seed(8)
NN <- 500
eta1 <- rnorm(NN); eta2 <- 0.6 * eta1 + rnorm(NN, 0, 0.8)
Yf <- data.frame(
  z1 = eta1 + rnorm(NN, 0, .5), z2 = 0.8 * eta1 + rnorm(NN, 0, .5),
  z3 = eta2 + rnorm(NN, 0, .5), z4 = 0.7 * eta2 + rnorm(NN, 0, .5)
)
Lf <- matrix(0, 4, 2); Lf[1:2, 1] <- 1; Lf[3:4, 2] <- 1
Bf <- matrix(0, 2, 2); Bf[2, 1] <- 1
m_fs <- suppressWarnings(
  runmodel(lvm(Yf, lambda = Lf, beta = Bf, identification = "variance"))
)
fs <- factorscores(Yf, m_fs, method = "regression")
L <- getmatrix(m_fs, "lambda"); Bm <- getmatrix(m_fs, "beta")
Psi <- getmatrix(m_fs, "sigma_zeta"); Th <- getmatrix(m_fs, "sigma_epsilon")
nu <- getmatrix(m_fs, "nu"); nu_eta <- getmatrix(m_fs, "nu_eta")
IminB <- solve(diag(2) - Bm)
SigEta <- IminB %*% Psi %*% t(IminB)
muEta <- as.vector(IminB %*% nu_eta)
SigY <- L %*% SigEta %*% t(L) + Th
muY <- as.vector(nu + L %*% muEta)
manual <- t(apply(Yf, 1, function(y) muEta + SigEta %*% t(L) %*% solve(SigY, y - muY)))
expect_true(max(abs(as.matrix(fs[, 1:2]) - manual)) < 1e-8)

## ---- latentgrowth: nonzero-mean LGC converges and matches lavaan ----
set.seed(7)
Ng <- 400
icept <- rnorm(Ng, 2, 1); slope <- rnorm(Ng, 0.5, 0.4)
lgc <- as.data.frame(sapply(0:3, function(tt) icept + slope * tt + rnorm(Ng, 0, 0.6)))
colnames(lgc) <- paste0("y", 1:4)
design_lgc <- matrix(paste0("y", 1:4), nrow = 1)
rownames(design_lgc) <- "y"
m_lg <- suppressWarnings(runmodel(latentgrowth(vars = design_lgc, data = lgc)))
expect_true(is.finite(m_lg@fitmeasures$logl))
expect_true(abs(m_lg@fitmeasures$logl - (-2220.18073)) < 0.01)
if (requireNamespace("lavaan", quietly = TRUE)) {
  lav_g <- lavaan::growth(
    "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4\ns =~ 0*y1 + 1*y2 + 2*y3 + 3*y4", lgc
  )
  lav_g_ll <- unname(lavaan::fitMeasures(lav_g)["logl"])
  expect_true(abs(m_lg@fitmeasures$logl - lav_g_ll) < 1e-3)
}

## ---- latentgrowth: nondefault time scores survive identification ----
# (0.15.29 fix: the identifier must not overwrite fixed nonzero loadings.)
m_lg2 <- suppressWarnings(latentgrowth(vars = design_lgc, data = lgc,
                                       time = c(0, 2, 4, 6)))
pt_lg <- m_lg2@parameters
slope_loadings <- pt_lg[pt_lg$matrix == "lambda" & pt_lg$col == 2, ]
expect_equal(sort(slope_loadings$est), c(0, 2, 4, 6))
expect_true(all(slope_loadings$fixed))

## ---- Ising: stray response values are rejected with a clear error ----
expect_error(
  Ising(data.frame(a = c(0, 1, 2, 0), b = c(0, 1, 0, 1)), responses = c(0L, 1L)),
  pattern = "not in 'responses'"
)

## ---- correlation input + mean structure is rejected for ML ----
expect_error(
  suppressWarnings(suppressMessages(
    ggm(covs = psych::Thurstone, means = rep(0, 9), nobs = 213)
  )),
  pattern = "corinput"
)

## ---- stepup() with ordinal data (GitHub issue #54) ----
# (0.16.8 fix: an EPC starting value overshooting into an improper region
# trapped the optimizer on the 1e20 penalty plateau; stepup() then crashed
# with "missing value where TRUE/FALSE needed" (NA BIC under DWLS) or
# silently rejected genuinely improving parameters.)
# The ordinal machinery calls the pbv package's C++ API; on some check setups
# (e.g. win-builder R-devel) the installed pbv binary lacks the Rcpp validation
# hook ("_pbv_RcppExport_validate not provided by package 'pbv'"), so any
# ordinal-data call errors. Probe the linkage first and skip if unusable:
pbv_usable <- tryCatch({
  psychonetrics:::covPrepare_cpp(
    cbind(c(1, 1, 2, 2, 1, 2, 2, 1), c(1, 2, 2, 1, 1, 2, 1, 2)),
    c(TRUE, TRUE), WLSweights = FALSE, verbose = FALSE)
  TRUE
}, error = function(e) FALSE)

if (pbv_usable) {
set.seed(1)
n_ord <- 500
eta_ord <- cbind(rnorm(n_ord), rnorm(n_ord))
eta_ord[, 2] <- 0.5 * eta_ord[, 1] + sqrt(1 - 0.25) * eta_ord[, 2]
lam_ord <- matrix(0, 8, 2)
lam_ord[1:4, 1] <- runif(4, .6, .9)
lam_ord[5:8, 2] <- runif(4, .6, .9)
err_ord <- matrix(rnorm(n_ord * 8), n_ord, 8) * sqrt(.5)
# Strong residual dependence between items 1 and 5:
shared <- rnorm(n_ord)
err_ord[, 1] <- err_ord[, 1] + sqrt(.3) * shared
err_ord[, 5] <- err_ord[, 5] + sqrt(.3) * shared
y_ord <- eta_ord %*% t(lam_ord) + err_ord
dat_ord <- as.data.frame(apply(y_ord, 2, function(x)
  as.numeric(cut(x, breaks = c(-Inf, quantile(x, c(.3, .6, .85)), Inf)))))
names(dat_ord) <- paste0("item", 1:8)
Lam_ord <- matrix(0, 8, 2)
Lam_ord[1:4, 1] <- 1
Lam_ord[5:8, 2] <- 1

m_ord <- suppressWarnings(suppressMessages(
  runmodel(rnm(dat_ord, lambda = Lam_ord, ordered = TRUE))
))
# stepup() must not error, and must recover the true residual edge:
m_ord_su <- suppressWarnings(suppressMessages(
  stepup(m_ord, criterion = "bic")
))
pt_ord <- m_ord_su@parameters
edge15 <- pt_ord[pt_ord$matrix == "omega_epsilon" &
                 pt_ord$row == 5 & pt_ord$col == 1, ]
expect_false(edge15$fixed)
# The refit from the EPC start must reach the proper optimum (not the
# improper-region penalty plateau):
expect_true(m_ord_su@objective < 1e10)
expect_true(m_ord_su@fitmeasures$chisq < m_ord@fitmeasures$chisq)
}

## ---- stepup(): EPC starting value outside the bounds of a parameter ----
# (0.18 fix: the expected parameter change of a partial correlation can be
# larger than 1. The implied covariance matrix is still positive definite at
# that value, the optimizer started on the bound, ended in a worse solution,
# and stepup() did not add the edge.)
# Population covariance matrix of two factors with four indicators each and
# one residual partial correlation (between y1 and y2):
Lam_epc <- matrix(0, 8, 2); Lam_epc[1:4, 1] <- 1; Lam_epc[5:8, 2] <- 1
Theta_epc <- diag(8); Theta_epc[2, 2] <- 1.36; Theta_epc[1, 2] <- Theta_epc[2, 1] <- 0.6
Sigma_epc <- 0.64 * Lam_epc %*% matrix(c(1, .4, .4, 1), 2, 2) %*% t(Lam_epc) + Theta_epc
dimnames(Sigma_epc) <- list(paste0("y", 1:8), paste0("y", 1:8))
m_epc <- suppressWarnings(suppressMessages(
  runmodel(rnm(covs = Sigma_epc, nobs = 500, covtype = "ML", lambda = Lam_epc))
))
edge21 <- function(m){
  m@parameters[m@parameters$matrix == "omega_epsilon" & m@parameters$row == 2 & m@parameters$col == 1, ]
}
# The expected value of the edge is not a partial correlation:
expect_true(edge21(m_epc)$est + edge21(m_epc)$epc > 1)
# freepar() then uses a small starting value:
expect_true(abs(edge21(freepar(m_epc, "omega_epsilon", 2, 1, verbose = FALSE))$est) < 0.01)
# A starting value that is given by hand is not changed:
expect_equal(edge21(freepar(m_epc, "omega_epsilon", 2, 1, start = 0.3, verbose = FALSE))$est, 0.3)
# stepup() adds the edge and recovers the population value:
m_epc_su <- suppressWarnings(suppressMessages(
  stepup(m_epc, matrices = "omega_epsilon")
))
expect_false(edge21(m_epc_su)$fixed)
expect_equal(m_epc_su@fitmeasures$df, m_epc@fitmeasures$df - 1)
expect_true(m_epc_su@fitmeasures$chisq < 1e-2)
Kappa_epc <- solve(Theta_epc)
expect_true(abs(edge21(m_epc_su)$est + Kappa_epc[1, 2] / sqrt(Kappa_epc[1, 1] * Kappa_epc[2, 2])) < 1e-2)

## ---- satMethodUsed / satLL_analytic tracking (GitHub issue #53) ----
# (0.16.8: the saturated-LL method ultimately used is recorded in
# @baseline_saturated$satMethodUsed and as the numeric fit measure
# satLL_analytic; saturated = "model" no longer silently falls back.)
bfi_53 <- psych::bfi[1:200, 1:5]
m_53 <- suppressWarnings(runmodel(ggm(bfi_53, estimator = "FIML")))
expect_equal(m_53@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_53@fitmeasures$satLL_analytic, 0)
m_53a <- suppressWarnings(
  runmodel(ggm(bfi_53, estimator = "FIML"), saturated = "analytic")
)
expect_equal(m_53a@baseline_saturated$satMethodUsed, "analytic")
expect_equal(m_53a@fitmeasures$satLL_analytic, 1)
# Force the fallback by sabotaging the numerically fitted saturated model:
sat_53 <- m_53@baseline_saturated$saturated
sel_off <- sat_53@parameters$matrix == "omega" &
  sat_53@parameters$row != sat_53@parameters$col
sel_dia <- sat_53@parameters$matrix == "delta" &
  sat_53@parameters$row == sat_53@parameters$col
sat_53@parameters$est[sel_off] <- 0
sat_53@parameters$est[sel_dia] <- 100
m_53f <- m_53
m_53f@baseline_saturated$saturated <- sat_53
m_53f@baseline_saturated$satMethod <- "default"
# (0.18: a saturated model that ends below the model is first re-estimated from
# the model-implied means and covariances, so this no longer needs the fallback.)
w_53f <- character(0)
m_53f <- withCallingHandlers(addfit(m_53f), warning = function(w){
  w_53f <<- c(w_53f, conditionMessage(w)); invokeRestart("muffleWarning")
})
expect_false(any(grepl("analytical saturated LL", w_53f)))
expect_equal(m_53f@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_53f@fitmeasures$satLL_analytic, 0)
expect_true(abs(m_53f@fitmeasures$chisq - m_53@fitmeasures$chisq) < 1e-3)
# The fallback remains for a saturated model that cannot reproduce the model
# (here the baseline model is put in its place):
m_53b <- m_53
m_53b@baseline_saturated$saturated <- m_53@baseline_saturated$baseline
m_53b@baseline_saturated$satMethod <- "default"
w_53b <- character(0)
m_53b <- withCallingHandlers(addfit(m_53b), warning = function(w){
  w_53b <<- c(w_53b, conditionMessage(w)); invokeRestart("muffleWarning")
})
expect_true(any(grepl("analytical saturated LL", w_53b)))
expect_equal(m_53b@baseline_saturated$satMethodUsed, "analytic_fallback")
expect_equal(m_53b@fitmeasures$satLL_analytic, 1)
# saturated = "model": documented as numeric with NO fallback:
m_53m <- m_53
m_53m@baseline_saturated$saturated <- sat_53
m_53m@baseline_saturated$satMethod <- "model"
m_53m <- suppressWarnings(addfit(m_53m))
expect_equal(m_53m@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_53m@fitmeasures$satLL_analytic, 0)

## ---- Saturated model re-estimated from the model-implied moments (0.18) ----
# A model with misfit of which the saturated model stopped far too early: the
# chi-square of the model must be recovered, without the analytic fallback.
data("NA2020", package = "psychonetrics")
onef_r <- function(d) suppressWarnings(suppressMessages(
  runmodel(lvm(d, lambda = matrix(1, ncol(d), 1), identification = "variance"))
))
m_re <- onef_r(NA2020)
sat_re <- m_re@baseline_saturated$saturated
sel_re <- sat_re@parameters$matrix != "mu" & sat_re@parameters$row != sat_re@parameters$col
sat_re@parameters$est[sel_re] <- 0
m_re2 <- m_re
m_re2@baseline_saturated$saturated <- sat_re
m_re2@baseline_saturated$satMethod <- "default"
w_re <- character(0)
m_re2 <- withCallingHandlers(addfit(m_re2), warning = function(w){
  w_re <<- c(w_re, conditionMessage(w)); invokeRestart("muffleWarning")
})
expect_false(any(grepl("analytical saturated LL", w_re)))
expect_equal(m_re2@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_re2@fitmeasures$chisq, m_re@fitmeasures$chisq, tolerance = 1e-4)
expect_equal(m_re2@fitmeasures$df, m_re@fitmeasures$df)

# A just-identified model with missing data has a chi-square of zero:
w_ji <- character(0)
m_ji <- withCallingHandlers(suppressMessages(
  runmodel(lvm(NA2020[, 1:3], lambda = matrix(1, 3, 1), identification = "variance"))
), warning = function(w){
  w_ji <<- c(w_ji, conditionMessage(w)); invokeRestart("muffleWarning")
})
expect_false(any(grepl("analytical saturated LL", w_ji)))
expect_equal(m_ji@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_ji@fitmeasures$df, 0)
expect_true(abs(m_ji@fitmeasures$chisq) < 1e-3)

## ---- FIML fit measures do not depend on the units of the data (0.18) ----
# Multiplying the data by a constant c lowers every log-likelihood by
# (number of observed values) * log(c) and leaves the chi-squares unchanged,
# also when the determinant of the variance-covariance matrix is large.
nobs_u <- sum(!is.na(NA2020))
m_u2 <- onef_r(NA2020 * 2)
expect_equal(m_u2@fitmeasures$logl, m_re@fitmeasures$logl - nobs_u * log(2), tolerance = 1e-6)
expect_equal(m_u2@fitmeasures$unrestricted.logl,
             m_re@fitmeasures$unrestricted.logl - nobs_u * log(2), tolerance = 1e-6)
expect_equal(m_u2@fitmeasures$baseline.logl,
             m_re@fitmeasures$baseline.logl - nobs_u * log(2), tolerance = 1e-6)
expect_equal(m_u2@fitmeasures$chisq, m_re@fitmeasures$chisq, tolerance = 1e-4)
expect_equal(m_u2@fitmeasures$baseline.chisq, m_re@fitmeasures$baseline.chisq, tolerance = 1e-4)
expect_equal(m_u2@fitmeasures$cfi, m_re@fitmeasures$cfi, tolerance = 1e-4)
expect_equal(m_u2@fitmeasures$bic, m_re@fitmeasures$bic + 2 * nobs_u * log(2), tolerance = 1e-6)
vc_u1 <- suppressWarnings(suppressMessages(runmodel(varcov(NA2020))))
vc_u10 <- suppressWarnings(suppressMessages(runmodel(varcov(NA2020 * 10))))
expect_equal(vc_u10@fitmeasures$logl, vc_u1@fitmeasures$logl - nobs_u * log(10), tolerance = 1e-5)

