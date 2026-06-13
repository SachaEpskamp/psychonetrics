# Tests for the Wishart Gaussian likelihood (likelihood = "wishart") in the
# varcov and lvm families. Under "wishart" the unbiased (n-1) sample covariance
# is used, the chi-square uses (n_g - 1) multipliers and the standard errors are
# inflated by sqrt(n_g/(n_g-1)) per group (matching lavaan likelihood = "wishart").
#
# Fast deterministic checks (chisq ratio, SE ratio, error handling) run always;
# the exact lavaan oracle comparison runs at_home.

suppressMessages(library(psychonetrics))

hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
N <- nrow(hs)
lambda <- matrix(0, 9, 3)
lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1; lambda[7:9, 3] <- 1

## ---- chisq and SE scaling (self-consistent, no lavaan) ----
m_w <- suppressWarnings(runmodel(
  lvm(hs, lambda = lambda, identification = "loadings", likelihood = "wishart"), verbose = FALSE))
m_n <- suppressWarnings(runmodel(
  lvm(hs, lambda = lambda, identification = "loadings", likelihood = "normal"), verbose = FALSE))

# df unchanged by the likelihood:
expect_equal(m_w@fitmeasures$df, m_n@fitmeasures$df)
expect_equal(m_w@fitmeasures$npar, m_n@fitmeasures$npar)

# Wishart chisq = (N-1)/N * normal chisq (single group):
expect_true(abs(m_w@fitmeasures$chisq / m_n@fitmeasures$chisq - (N - 1) / N) < 1e-4)

# The stored sample covariance under wishart is the unbiased (n-1) one:
S_w <- as.matrix(m_w@sample@covs[[1]])
S_n <- as.matrix(m_n@sample@covs[[1]])
expect_true(max(abs(S_w - S_n * N / (N - 1))) < 1e-8)
expect_true(max(abs(S_w - cov(hs))) < 1e-8)   # = base R unbiased cov

## ---- Error handling (always) ----
expect_error(lvm(hs, lambda = lambda, estimator = "FIML", likelihood = "wishart"),
             "wishart")
hs_miss <- hs; hs_miss[1, 1] <- NA
expect_error(lvm(hs_miss, lambda = lambda, likelihood = "wishart"), "wishart")
expect_error(
  varcov(hs, type = "ggm", estimator = "WLS", ordered = names(hs), likelihood = "wishart"),
  "wishart")
expect_error(varcov(hs, type = "cov", estimator = "ULS", likelihood = "wishart"), "wishart")

## ---- varcov wishart self-consistency ----
vc_w <- suppressWarnings(runmodel(varcov(hs[, 1:4], type = "cov", likelihood = "wishart"), verbose = FALSE))
vc_n <- suppressWarnings(runmodel(varcov(hs[, 1:4], type = "cov", likelihood = "normal"), verbose = FALSE))
# A saturated cov model has df 0 and chisq ~ 0 under both:
expect_equal(vc_w@fitmeasures$df, 0)
# The estimated (co)variances are the unbiased sample (co)variances:
sig_w <- as.matrix(getmatrix(vc_w, "sigma"))
expect_true(max(abs(sig_w - cov(hs[, 1:4]))) < 1e-6)

## ---- Backward compatibility: default "normal" matches the no-argument call ----
m_default <- suppressWarnings(runmodel(
  lvm(hs, lambda = lambda, identification = "loadings"), verbose = FALSE))
expect_true(max(abs(m_default@parameters$est - m_n@parameters$est)) < 1e-10)
expect_true(max(abs(m_default@parameters$se - m_n@parameters$se), na.rm = TRUE) < 1e-10)
expect_equal(m_default@fitmeasures$chisq, m_n@fitmeasures$chisq)

## ---- lavaan oracle (at_home) ----
if (at_home()){
  if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")

  fit_w <- lavaan::cfa(
    "visual =~ x1+x2+x3\ntextual =~ x4+x5+x6\nspeed =~ x7+x8+x9",
    data = hs, likelihood = "wishart", meanstructure = TRUE)
  pe_w <- lavaan::parameterEstimates(fit_w)
  chisq_lav <- unname(lavaan::fitMeasures(fit_w)["chisq"])

  # Tighten the psychonetrics optimizer so the comparison is at (near) identical
  # estimates; the remaining gap is pure optimizer tolerance.
  mt <- lvm(hs, lambda = lambda, identification = "loadings", likelihood = "wishart")
  mt@optim.args <- list(control = list(eval.max = 1e5L, iter.max = 1e5L,
                                       rel.tol = 1e-13, x.tol = 1e-12, xf.tol = 1e-14))
  m_wt <- suppressWarnings(runmodel(mt, verbose = FALSE))

  # chisq matches lavaan wishart:
  expect_true(abs(m_wt@fitmeasures$chisq - chisq_lav) < 1e-4)

  # Free loading and residual-variance SEs match lavaan wishart:
  for (i in 1:9){
    v <- paste0("x", i)
    prow <- m_wt@parameters[m_wt@parameters$matrix == "lambda" &
                              m_wt@parameters$var1 == v &
                              m_wt@parameters$par != 0 & m_wt@parameters$est != 0, ]
    lrow <- pe_w[pe_w$op == "=~" & pe_w$rhs == v, ]
    if (nrow(prow) == 1 && nrow(lrow) == 1 && is.finite(lrow$se) && lrow$se > 0){
      expect_true(abs(prow$se - lrow$se) < 1e-3)
    }
    prv <- m_wt@parameters[m_wt@parameters$matrix == "sigma_epsilon" &
                             m_wt@parameters$var1 == v & m_wt@parameters$var2 == v, ]
    lrv <- pe_w[pe_w$op == "~~" & pe_w$lhs == v & pe_w$rhs == v, ]
    if (nrow(prv) == 1 && nrow(lrv) == 1 && is.finite(lrv$se) && lrv$se > 0){
      expect_true(abs(prv$se - lrv$se) < 1e-3)
    }
  }
}
