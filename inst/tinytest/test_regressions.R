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
expect_warning(m_53f <- addfit(m_53f), pattern = "analytical saturated LL")
expect_equal(m_53f@baseline_saturated$satMethodUsed, "analytic_fallback")
expect_equal(m_53f@fitmeasures$satLL_analytic, 1)
# saturated = "model": documented as numeric with NO fallback:
m_53m <- m_53
m_53m@baseline_saturated$saturated <- sat_53
m_53m@baseline_saturated$satMethod <- "model"
m_53m <- suppressWarnings(addfit(m_53m))
expect_equal(m_53m@baseline_saturated$satMethodUsed, "numeric")
expect_equal(m_53m@fitmeasures$satLL_analytic, 0)
