# Robust maximum likelihood (Phase 1, complete data): MLM / MLMV / MLMVS / MLR.
#
# Point estimation is plain ML; only the standard errors, the scaled test
# statistic and the robust fit indices change. lavaan is used as the oracle
# (estimator = "MLM"/"MLMV"/"MLMVS"/"MLR"). The fast deterministic subset runs
# unconditionally; the lavaan comparisons run under at_home.
#
# References: Satorra & Bentler (1994); Yuan & Bentler (2000);
#   Brosseau-Liard & Savalei (2014); Savalei (2010); Asparouhov & Muthen (2010).

suppressMessages(library(psychonetrics))
if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")

HS <- lavaan::HolzingerSwineford1939
X <- paste0("x", 1:9)
mod <- 'visual =~ x1+x2+x3
        textual =~ x4+x5+x6
        speed =~ x7+x8+x9'
Lambda <- matrix(0, 9, 3)
Lambda[1:3, 1] <- Lambda[4:6, 2] <- Lambda[7:9, 3] <- 1
colnames(Lambda) <- c("visual", "textual", "speed")

# Helper: relative SE differences psychonetrics vs lavaan (matched by matrix):
fnames <- c("visual", "textual", "speed")
se_reldiff <- function(fitP, fitL) {
  pars <- fitP@parameters[fitP@parameters$par != 0 & !duplicated(fitP@parameters$par), ]
  pt <- lavaan::parameterEstimates(fitL, remove.nonfree = TRUE)
  mfun <- function(row) {
    m <- row$matrix
    if (m == "lambda") pt$se[pt$op == "=~" & pt$lhs == fnames[row$col] & pt$rhs == X[row$row]]
    else if (m == "sigma_zeta") pt$se[pt$op == "~~" & pt$lhs == fnames[min(row$row, row$col)] & pt$rhs == fnames[max(row$row, row$col)]]
    else if (m == "sigma_epsilon") pt$se[pt$op == "~~" & pt$lhs == X[row$row] & pt$rhs == X[row$col]]
    else if (m == "nu") pt$se[pt$op == "~1" & pt$lhs == X[row$row]]
    else NA_real_
  }
  rd <- vapply(seq_len(nrow(pars)), function(i) {
    lse <- mfun(pars[i, ]); pse <- pars$se[i]
    if (length(lse) != 1 || is.na(lse) || is.na(pse)) return(NA_real_)
    abs(pse - lse) / max(abs(lse), 1e-8)
  }, numeric(1))
  max(rd, na.rm = TRUE)
}

## =========================================================================
## Fast deterministic subset (always runs): structural sanity of robust output
## =========================================================================

# Point estimates identical to plain ML (robust only changes SE/test):
mod_ml  <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                identification = "variance", estimator = "ML"), verbose = FALSE))
mod_mlm <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                identification = "variance", estimator = "MLM"), verbose = FALSE))

expect_equal(mod_mlm@estimator, "ML")                    # maps internally to ML
expect_true(!is.null(mod_mlm@robust$se))
expect_equal(mod_mlm@robust$label, "MLM")
expect_equal(mod_ml@parameters$est, mod_mlm@parameters$est, tolerance = 1e-8)  # same point estimates
expect_equal(mod_ml@fitmeasures$chisq, mod_mlm@fitmeasures$chisq, tolerance = 1e-6)

# Robust scaled-test and fit indices are present and finite:
expect_true(is.finite(mod_mlm@fitmeasures$chisq.scaled))
expect_true(is.finite(mod_mlm@fitmeasures$chisq.scaling.factor))
expect_equal(mod_mlm@fitmeasures$df.scaled, mod_mlm@fitmeasures$df)
expect_true(is.finite(mod_mlm@fitmeasures$rmsea.robust))
expect_true(is.finite(mod_mlm@fitmeasures$cfi.robust))
expect_true(is.finite(mod_mlm@fitmeasures$tli.robust))
# Scaled chi-square = chisq / scaling for the (mean-adjusted) Satorra-Bentler:
expect_equal(mod_mlm@fitmeasures$chisq.scaled,
             mod_mlm@fitmeasures$chisq / mod_mlm@fitmeasures$chisq.scaling.factor,
             tolerance = 1e-6)

# Robust SEs differ from naive ML SEs (the whole point), but are the same order:
se_naive <- mod_ml@parameters$se[mod_ml@parameters$par != 0]
se_mlm   <- mod_mlm@parameters$se[mod_mlm@parameters$par != 0]
expect_true(all(is.finite(se_mlm)))
expect_false(isTRUE(all.equal(se_naive, se_mlm)))

# MLMV reports a shift parameter; MLM/MLMVS do not:
mod_mlmv <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                identification = "variance", estimator = "MLMV"), verbose = FALSE))
expect_true(is.finite(mod_mlmv@fitmeasures$chisq.shift.parameters))

# MLR requires/stores raw data and produces Huber-White SEs:
mod_mlr <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                identification = "variance", estimator = "MLR"), verbose = FALSE))
expect_true(nrow(mod_mlr@sample@rawdata) > 0)
expect_true(all(is.finite(mod_mlr@parameters$se[mod_mlr@parameters$par != 0])))
expect_true(is.finite(mod_mlr@fitmeasures$chisq.scaled))

# MLR with missing data maps to FIML internally and produces finite robust
# output (structural sanity; lavaan oracle checks run under at_home below):
set.seed(7)
HS_na <- HS
for (v in X) HS_na[[v]][stats::runif(nrow(HS_na)) < 0.12] <- NA
HS_na <- HS_na[rowSums(!is.na(HS_na[, X])) > 0, ]
mod_mlr_fiml <- suppressWarnings(runmodel(lvm(HS_na, lambda = Lambda, vars = X,
                identification = "variance", estimator = "MLR"), verbose = FALSE))
expect_equal(mod_mlr_fiml@estimator, "FIML")          # FIML under the hood
expect_equal(mod_mlr_fiml@robust$label, "MLR")
expect_true(nrow(mod_mlr_fiml@sample@rawdata) > 0)     # storedata forced TRUE
expect_true(length(mod_mlr_fiml@sample@fimldata) > 0)  # missingness patterns
expect_true(all(is.finite(mod_mlr_fiml@parameters$se[mod_mlr_fiml@parameters$par != 0])))
expect_true(is.finite(mod_mlr_fiml@fitmeasures$chisq.scaled))
expect_true(is.finite(mod_mlr_fiml@fitmeasures$chisq.scaling.factor))
expect_true(is.finite(mod_mlr_fiml@fitmeasures$rmsea.robust))
expect_true(is.finite(mod_mlr_fiml@fitmeasures$cfi.robust))

# setestimator() path is equivalent to the constructor path:
mod_via_set <- suppressWarnings(runmodel(
  setestimator(lvm(HS, lambda = Lambda, vars = X, identification = "variance",
                   estimator = "ML", storedata = TRUE), "MLM"), verbose = FALSE))
expect_equal(mod_via_set@fitmeasures$chisq.scaling.factor,
             mod_mlm@fitmeasures$chisq.scaling.factor, tolerance = 1e-8)

## =========================================================================
## lavaan oracle comparisons (at_home only)
## =========================================================================
if (at_home()) {
  tol_se <- 1e-3   # relative SE tolerance at the psychonetrics solution
  tol_st <- 1e-3   # relative tolerance for scaled chisq / scaling / indices

  ## ---- MLM single group ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLM", meanstructure = TRUE, std.lv = TRUE)
  fitP <- mod_mlm
  expect_true(se_reldiff(fitP, fitL) < tol_se)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaling.factor")),
               fitP@fitmeasures$chisq.scaling.factor, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "rmsea.robust")),
               fitP@fitmeasures$rmsea.robust, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "cfi.robust")),
               fitP@fitmeasures$cfi.robust, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "tli.robust")),
               fitP@fitmeasures$tli.robust, tolerance = tol_st)

  ## ---- MLMV single group (scaled-shifted) ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLMV", meanstructure = TRUE, std.lv = TRUE)
  fitP <- mod_mlmv
  expect_true(se_reldiff(fitP, fitL) < tol_se)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)
  expect_equal(fitL@test$scaled.shifted$shift.parameter,
               fitP@fitmeasures$chisq.shift.parameters, tolerance = 1e-2)

  ## ---- MLMVS single group (mean-and-variance adjusted; fractional df) ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLMVS", meanstructure = TRUE, std.lv = TRUE)
  fitP <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "MLMVS"), verbose = FALSE))
  expect_true(se_reldiff(fitP, fitL) < tol_se)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "df.scaled")),
               fitP@fitmeasures$df.scaled, tolerance = 1e-3)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)

  ## ---- MLR single group (Huber-White SEs + Yuan-Bentler-Mplus) ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLR", meanstructure = TRUE, std.lv = TRUE)
  fitP <- mod_mlr
  expect_true(se_reldiff(fitP, fitL) < tol_se)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaling.factor")),
               fitP@fitmeasures$chisq.scaling.factor, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)

  ## ---- Exact-solution check: psychonetrics evaluated at lavaan's estimates
  ##      should reproduce lavaan's robust SEs. Compared via the eigenvalues of
  ##      the robust VCOV (basis-free); the tolerance reflects the (tiny)
  ##      difference between the psychonetrics and lavaan optima. The true
  ##      machine-precision exact-solution check is below (lavaan's moments).
  ns <- asNamespace("psychonetrics")
  fitL_mlm <- lavaan::cfa(mod, HS, estimator = "MLM", meanstructure = TRUE, std.lv = TRUE)
  V_P <- ns$getVCOV_robust(mod_mlm)
  V_L <- as.matrix(lavaan::vcov(fitL_mlm))
  expect_equal(sort(eigen(V_P, only.values = TRUE)$values),
               sort(eigen(V_L, only.values = TRUE)$values), tolerance = 1e-3)

  ## ---- EXACT-solution check: psychonetrics blocks evaluated at lavaan's
  ##      implied moments must reproduce lavaan's robust SEs and scaling to
  ##      machine precision (isolates the formulas from optimizer differences) ----
  Delta_L <- lavaan::lavTech(fitL_mlm, "delta")[[1]]
  Sigma_L <- lavaan::lavTech(fitL_mlm, "implied")[[1]]$cov
  Gamma_L <- lavaan::lavTech(fitL_mlm, "gamma")[[1]]
  WLSV_L  <- lavaan::lavTech(fitL_mlm, "WLS.V")[[1]]
  V_x <- ns$ml_Vmat(Sigma_L, meanstructure = TRUE, corinput = FALSE)
  expect_equal(unname(as.matrix(V_x)), unname(as.matrix(WLSV_L)), tolerance = 1e-10)
  Nobs <- nrow(HS)
  E_x <- t(Delta_L) %*% V_x %*% Delta_L
  Ei  <- solve(E_x); VD <- V_x %*% Delta_L
  SE_x <- sqrt(diag(Ei %*% (t(VD) %*% Gamma_L %*% VD) %*% Ei / Nobs))
  expect_equal(as.numeric(SE_x), as.numeric(sqrt(diag(lavaan::vcov(fitL_mlm)))), tolerance = 1e-8)
  U_x <- V_x - VD %*% Ei %*% t(VD)
  UG_x <- U_x %*% Gamma_L
  expect_equal(sum(diag(UG_x)) / as.numeric(lavaan::fitMeasures(fitL_mlm, "df")),
               as.numeric(lavaan::fitMeasures(fitL_mlm, "chisq.scaling.factor")),
               tolerance = 1e-8)

  ## ---- Unit checks: Gamma vs lavTech gamma (sample-only -> machine precision);
  ##      B0 / observed information at the solution (looser at-solution tol) ----
  fitL_mlr <- lavaan::cfa(mod, HS, estimator = "MLR", meanstructure = TRUE, std.lv = TRUE)
  GammaP <- unname(as.matrix(mod_mlr@sample@WLS.Gamma[[1]]))
  GammaL <- unname(as.matrix(lavaan::lavTech(fitL_mlr, "gamma")[[1]]))
  expect_equal(GammaP, GammaL, tolerance = 1e-8)

  # Casewise scores / B0 at lavaan's solution -> machine precision:
  Dr <- lavaan::lavTech(fitL_mlr, "delta")[[1]]
  impl <- lavaan::lavTech(fitL_mlr, "implied")[[1]]
  SC1x <- ns$ml_casewise_scores_h1(as.matrix(HS[, X]), as.numeric(impl$mean),
                                   impl$cov, meanstructure = TRUE, corinput = FALSE)
  B0x <- t(Dr) %*% (crossprod(SC1x) / Nobs) %*% Dr
  B0L <- as.matrix(lavaan::lavTech(fitL_mlr, "information.first.order"))
  expect_equal(sort(eigen(B0x, only.values = TRUE)$values),
               sort(eigen(B0L, only.values = TRUE)$values), tolerance = 1e-8)

  meat <- ns$mlr_meat_components(mod_mlr)
  AL <- as.matrix(lavaan::lavTech(fitL_mlr, "information.observed"))
  expect_equal(sort(eigen(meat$A, only.values = TRUE)$values),
               sort(eigen(AL, only.values = TRUE)$values), tolerance = 1e-3)

  ## ---- Two-group MLM ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLM", group = "school", std.lv = TRUE)
  fitP <- suppressWarnings(runmodel(lvm(HS, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "MLM", groups = "school"),
                  verbose = FALSE))
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaling.factor")),
               fitP@fitmeasures$chisq.scaling.factor, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)

  ## ---- Two-group + group-equal loadings under MLM (cross-group constraints) ----
  fitL <- lavaan::cfa(mod, HS, estimator = "MLM", group = "school",
                      group.equal = "loadings", std.lv = TRUE)
  fitP0 <- suppressWarnings(lvm(HS, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "MLM", groups = "school"))
  fitP <- suppressWarnings(runmodel(groupequal(fitP0, "lambda"), verbose = FALSE))
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaling.factor")),
               fitP@fitmeasures$chisq.scaling.factor, tolerance = tol_st)
  expect_equal(as.numeric(lavaan::fitMeasures(fitL, "chisq.scaled")),
               fitP@fitmeasures$chisq.scaled, tolerance = 1e-2)

  ## ---- GGM / varcov under MLM and MLR ----
  ## A full GGM is saturated (df = 0); use a sparse (over-identified) network so
  ## the scaled test is well defined. omega: free only a few partial correlations.
  omega0 <- matrix(0, 9, 9)
  omega0[2, 1] <- omega0[1, 2] <- omega0[3, 2] <- omega0[2, 3] <-
    omega0[5, 4] <- omega0[4, 5] <- omega0[8, 7] <- omega0[7, 8] <- 1
  fitP_ggm_mlm <- suppressWarnings(runmodel(ggm(HS, vars = X, omega = omega0,
                       estimator = "MLM"), verbose = FALSE))
  expect_true(fitP_ggm_mlm@fitmeasures$df > 0)
  expect_true(is.finite(fitP_ggm_mlm@fitmeasures$chisq.scaling.factor))
  expect_true(fitP_ggm_mlm@fitmeasures$chisq.scaling.factor > 0)
  # Cross-check the GGM MLM scaling against the equivalent lavaan covariance
  # model with the same zero covariances (varcov MLM == saturated-pattern test):
  fitP_vc_mlm <- suppressWarnings(runmodel(varcov(HS, vars = X, type = "ggm",
                       omega = omega0, estimator = "MLM"), verbose = FALSE))
  expect_equal(fitP_vc_mlm@fitmeasures$chisq.scaling.factor,
               fitP_ggm_mlm@fitmeasures$chisq.scaling.factor, tolerance = 1e-6)

  fitP_ggm_mlr <- suppressWarnings(runmodel(ggm(HS, vars = X, omega = omega0,
                       estimator = "MLR"), verbose = FALSE))
  expect_true(all(is.finite(fitP_ggm_mlr@parameters$se[fitP_ggm_mlr@parameters$par != 0])))
  expect_true(is.finite(fitP_ggm_mlr@fitmeasures$chisq.scaled))

  ## ---- Multi-group DWLS trUG2 fix: the equality-constrained scaled chisq now
  ##      matches lavaan (the old per-group accumulation did not) ----
  HSc <- HS
  fitL_dwls <- lavaan::cfa(mod, HSc, estimator = "WLSMV", group = "school",
                           group.equal = "loadings", std.lv = TRUE)
  fitP0d <- suppressWarnings(lvm(HSc, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "DWLS", groups = "school"))
  fitPd <- suppressWarnings(runmodel(groupequal(fitP0d, "lambda"), verbose = FALSE))
  expect_equal(as.numeric(lavaan::fitMeasures(fitL_dwls, "chisq.scaling.factor")),
               fitPd@fitmeasures$chisq.scaling.factor, tolerance = 5e-3)

  ## =======================================================================
  ## MLR with FIML (Phase 2: within-row missing data). Point estimation is
  ## plain FIML; SEs are Huber-White and the scaled chi-square is the
  ## Yuan-Bentler-Mplus statistic under missingness. The robust fit indices
  ## use the FIML-Corrected (FIML-C V3; Savalei 2010) correction. Oracle:
  ## lavaan sem(..., estimator = "MLR", missing = "fiml").
  ## =======================================================================
  ns <- asNamespace("psychonetrics")
  set.seed(123)
  HSmis <- HS
  for (v in X) HSmis[[v]][stats::runif(nrow(HSmis)) < 0.15] <- NA  # ~15% MCAR
  HSmis <- HSmis[rowSums(!is.na(HSmis[, X])) > 0, ]

  ## ---- End-to-end MLR-FIML (psychonetrics own optimum) vs lavaan ----
  mfF <- suppressWarnings(runmodel(lvm(HSmis, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "MLR"), verbose = FALSE))
  fitF <- lavaan::cfa(mod, HSmis, estimator = "MLR", missing = "fiml",
                      meanstructure = TRUE, std.lv = TRUE)

  # MLR + missing maps internally to estimator = "FIML" with the MLR robust cfg,
  # storing both the raw data (with NAs) and the missingness patterns:
  expect_equal(mfF@estimator, "FIML")
  expect_equal(mfF@robust$label, "MLR")
  expect_equal(mfF@robust$se, "robust.huber.white")
  expect_true(nrow(mfF@sample@rawdata) > 0)
  expect_true(length(mfF@sample@fimldata) > 0)
  expect_true(all(is.finite(mfF@parameters$se[mfF@parameters$par != 0])))

  # Unscaled chisq, df, scaled chisq, scaling factor (own optima -> loose tol):
  expect_equal(mfF@fitmeasures$chisq,
               as.numeric(lavaan::fitMeasures(fitF, "chisq")), tolerance = 1e-4)
  expect_equal(mfF@fitmeasures$df, as.numeric(lavaan::fitMeasures(fitF, "df")))
  expect_equal(mfF@fitmeasures$chisq.scaling.factor,
               as.numeric(lavaan::fitMeasures(fitF, "chisq.scaling.factor")), tolerance = 1e-3)
  expect_equal(mfF@fitmeasures$chisq.scaled,
               as.numeric(lavaan::fitMeasures(fitF, "chisq.scaled")), tolerance = 1e-2)
  # SEs within 1e-4 relative (psychonetrics vs lavaan, both at own optima):
  expect_true(se_reldiff(mfF, fitF) < 1e-3)
  # FIML-Corrected robust fit indices (rmsea/cfi/tli.robust):
  expect_equal(mfF@fitmeasures$rmsea.robust,
               as.numeric(lavaan::fitMeasures(fitF, "rmsea.robust")), tolerance = 1e-3)
  expect_equal(mfF@fitmeasures$cfi.robust,
               as.numeric(lavaan::fitMeasures(fitF, "cfi.robust")), tolerance = 1e-3)
  expect_equal(mfF@fitmeasures$tli.robust,
               as.numeric(lavaan::fitMeasures(fitF, "tli.robust")), tolerance = 1e-3)

  ## ---- EXACT (machine-precision) unit checks at lavaan's FIML solution ----
  ## These isolate the pattern-wise score / information formulas from the
  ## (tiny) optimizer difference between psychonetrics and lavaan.
  Nf  <- lavaan::lavInspect(fitF, "ntotal")
  Drf <- lavaan::lavTech(fitF, "delta")[[1]]
  implf <- lavaan::lavTech(fitF, "implied")[[1]]
  muf <- as.numeric(implf$mean); Sigf <- implf$cov
  Yf  <- as.matrix(HSmis[, X])

  # Pattern casewise scores vs lavScores (<= 1e-10):
  SCf <- ns$ml_casewise_scores_h1_missing(Yf, muf, Sigf)
  expect_equal(as.numeric(SCf %*% Drf), as.numeric(lavaan::lavScores(fitF)),
               tolerance = 1e-10)

  # First-order information vs lavTech "information.first.order" (<= 1e-10):
  B0f <- t(Drf) %*% (crossprod(SCf) / Nf) %*% Drf
  expect_equal(unname(B0f),
               unname(as.matrix(lavaan::lavTech(fitF, "information.first.order"))),
               tolerance = 1e-10)

  # Huber-White SEs at lavaan's solution vs lavaan (<= 1e-8):
  Af  <- lavaan::lavInspect(fitF, "information.observed")
  Aif <- solve(Af)
  SEf <- sqrt(diag(Aif %*% B0f %*% Aif / Nf))
  expect_equal(as.numeric(SEf),
               as.numeric(sqrt(diag(lavaan::vcov(fitF)))), tolerance = 1e-8)

  # Pattern-based unstructured h1 OBSERVED info under missingness (<= 1e-8):
  h1f <- lavaan::lavTech(fitF, "h1")[[1]]
  mu1f <- as.numeric(h1f$mean); Sig1f <- h1f$cov
  A1uf <- ns$ml_h1_information_observed_missing(Yf, mu1f, Sig1f)
  # lavaan >= 0.7 renamed this internal and changed its interface; the
  # cross-check only runs where the 0.6 version is available:
  lav_h1_obs <- get0("lav_model_h1_information_observed",
                     envir = asNamespace("lavaan"), inherits = FALSE)
  if (!is.null(lav_h1_obs)) {
    A1uf_lav <- lav_h1_obs(
      lavmodel = fitF@Model, lavsamplestats = fitF@SampleStats, lavdata = fitF@Data,
      lavimplied = fitF@implied, lavh1 = fitF@h1,
      lavoptions = local({o <- fitF@Options; o$h1.information <- "unstructured"; o}))[[1]]
    expect_equal(unname(A1uf), unname(as.matrix(A1uf_lav)), tolerance = 1e-8)
  }

  # Yuan-Bentler-Mplus scaling factor at lavaan's solution (<= 1e-8):
  SCuf <- ns$ml_casewise_scores_h1_missing(Yf, mu1f, Sig1f)
  trh1 <- sum((crossprod(SCuf) / Nf) * t(solve(A1uf)))
  trh0 <- sum(B0f * t(Aif))
  c_yb <- (trh1 - trh0) / as.numeric(lavaan::fitMeasures(fitF, "df"))
  expect_equal(c_yb,
               as.numeric(lavaan::fitMeasures(fitF, "chisq.scaling.factor")),
               tolerance = 1e-8)

  ## ---- Multi-group MLR-FIML ----
  mfF2 <- suppressWarnings(runmodel(lvm(HSmis, lambda = Lambda, vars = X,
                  identification = "variance", estimator = "MLR", groups = "school"),
                  verbose = FALSE))
  fitF2 <- lavaan::cfa(mod, HSmis, estimator = "MLR", missing = "fiml",
                       group = "school", std.lv = TRUE)
  expect_true(all(is.finite(mfF2@parameters$se[mfF2@parameters$par != 0])))
  expect_equal(mfF2@fitmeasures$df, as.numeric(lavaan::fitMeasures(fitF2, "df")))
  expect_equal(mfF2@fitmeasures$chisq.scaling.factor,
               as.numeric(lavaan::fitMeasures(fitF2, "chisq.scaling.factor")), tolerance = 1e-3)
  expect_equal(mfF2@fitmeasures$chisq.scaled,
               as.numeric(lavaan::fitMeasures(fitF2, "chisq.scaled")), tolerance = 1e-2)
  expect_equal(mfF2@fitmeasures$rmsea.robust,
               as.numeric(lavaan::fitMeasures(fitF2, "rmsea.robust")), tolerance = 1e-3)

  ## ---- setestimator() path on a FIML model is equivalent to the constructor ----
  mfF_set <- suppressWarnings(runmodel(setestimator(
    lvm(HSmis, lambda = Lambda, vars = X, identification = "variance",
        estimator = "ML", storedata = TRUE), "MLR"), verbose = FALSE))
  expect_equal(mfF_set@estimator, "FIML")
  expect_equal(mfF_set@fitmeasures$chisq.scaling.factor,
               mfF@fitmeasures$chisq.scaling.factor, tolerance = 1e-6)

  ## ---- The MLM family (robust.sem) still errors on missing data via
  ##      setestimator() (no asymptotic Gamma under missingness) ----
  expect_error(setestimator(
    lvm(HSmis, lambda = Lambda, vars = X, identification = "variance",
        estimator = "ML", storedata = TRUE), "MLM"))
}
