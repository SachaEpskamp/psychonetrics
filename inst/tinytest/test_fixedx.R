# Tests for fixed.x exogenous covariates (fixed_x = ...) in the varcov and lvm
# families. The means and mutual (co)variances of the declared exogenous
# variables are fixed to their sample values and excluded from npar / df,
# conditioning the model on x (matching lavaan's sem(..., fixed.x = TRUE);
# Rosseel, 2012). The reported log-likelihood is the conditional log-likelihood.
#
# Deterministic structural checks run always; the exact lavaan oracle
# comparisons (estimates, SEs, chisq, df, logl) run at_home.

suppressMessages(library(psychonetrics))

## ---- varcov fixed_x: structural checks (no lavaan) ----
dat <- lavaan::PoliticalDemocracy[, c("y1", "x1", "x2", "x3")]
m_fx <- suppressWarnings(runmodel(
  varcov(dat, type = "cov", vars = c("y1", "x1", "x2", "x3"),
         fixed_x = c("x1", "x2", "x3")), verbose = FALSE))
m_plain <- suppressWarnings(runmodel(varcov(dat, type = "cov"), verbose = FALSE))

# The x-block (3 means + 6 x-x (co)variances = 9 statistics) leaves npar:
expect_equal(m_fx@fitmeasures$npar, m_plain@fitmeasures$npar - 9)
# Saturated either way -> df 0:
expect_equal(m_fx@fitmeasures$df, 0)
# The fixed x-block parameters equal the sample moments:
pt <- m_fx@parameters
S <- as.matrix(m_fx@sample@covs[[1]])
mn <- as.numeric(m_fx@sample@means[[1]])
xidx <- 2:4
for (i in which(pt$matrix == "sigma" & pt$row %in% xidx & pt$col %in% xidx)){
  expect_true(pt$fixed[i]); expect_equal(pt$par[i], 0)
  expect_equal(pt$est[i], S[pt$row[i], pt$col[i]], tolerance = 1e-8)
}
for (i in which(pt$matrix == "mu" & pt$row %in% xidx)){
  expect_true(pt$fixed[i]); expect_equal(pt$est[i], mn[pt$row[i]], tolerance = 1e-8)
}
# The implied regression of y1 on x equals the OLS / SEM coefficients:
Sig <- as.matrix(getmatrix(m_fx, "sigma"))
b <- solve(Sig[xidx, xidx], Sig[xidx, 1])
ols <- coef(lm(y1 ~ x1 + x2 + x3, data = dat))[-1]
expect_equal(unname(b), unname(ols), tolerance = 1e-6)

## ---- Default fixed_x = character(0) is unchanged ----
m_def <- suppressWarnings(runmodel(varcov(dat, type = "cov"), verbose = FALSE))
expect_equal(m_def@fitmeasures$npar, 14L)   # 4 vars: 4 means + 10 (co)vars
expect_true(is.null(m_def@types$fixed_x) || length(m_def@types$fixed_x) == 0)

## ---- fixed_x error handling ----
expect_error(varcov(dat, type = "ggm", fixed_x = c("x1")), "type = 'cov'")
expect_error(
  varcov(dat, type = "cov", vars = c("y1","x1","x2","x3"), fixed_x = c("nope")),
  "not found")
expect_error(
  varcov(dat, type = "cov", vars = c("y1","x1","x2","x3"),
         fixed_x = c("y1","x1","x2","x3")),
  "at least one endogenous")

## ---- lvm fixed_x: single-indicator-latent SEM structural checks ----
vars <- c("y1", "y2", "y3", "y4", "x1", "x2", "x3")
dat2 <- lavaan::PoliticalDemocracy[, vars]
lambda <- matrix(0, 7, 4)
lambda[1:4, 1] <- 1; lambda[5, 2] <- 1; lambda[6, 3] <- 1; lambda[7, 4] <- 1
beta <- matrix(0, 4, 4); beta[1, 2] <- 1; beta[1, 3] <- 1; beta[1, 4] <- 1
se <- matrix(0, 7, 7); diag(se)[1:4] <- 1
lats <- c("dem60", "Lx1", "Lx2", "Lx3")
m_lvm_fx <- suppressWarnings(runmodel(
  lvm(dat2, lambda = lambda, beta = beta, sigma_epsilon = se, vars = vars,
      identification = "loadings", latents = lats, fixed_x = c("Lx1", "Lx2", "Lx3")),
  verbose = FALSE))
m_lvm_plain <- suppressWarnings(runmodel(
  lvm(dat2, lambda = lambda, beta = beta, sigma_epsilon = se, vars = vars,
      identification = "loadings", latents = lats), verbose = FALSE))
# x-block (3 latent vars + 3 latent covs + 3 observed x means = 9) leaves npar:
expect_equal(m_lvm_fx@fitmeasures$npar, m_lvm_plain@fitmeasures$npar - 9)
# df increases by 9 (statistics fixed) ... actually nobs and npar both drop by 9,
# so df is unchanged relative to plain (both condition on the same moments):
expect_equal(m_lvm_fx@fitmeasures$df, m_lvm_plain@fitmeasures$df)

# fixed_x must name latents that are single-indicator:
expect_error(
  lvm(dat2, lambda = lambda, beta = beta, sigma_epsilon = se, vars = vars,
      latents = lats, fixed_x = c("x1")),
  "latent")

## ---- lavaan oracle (at_home) ----
if (at_home()){
  if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")

  # (1) Saturated regression: estimates / df / npar / logl exact.
  fitL <- lavaan::sem("y1 ~ x1 + x2 + x3", data = dat, fixed.x = TRUE,
                      meanstructure = TRUE)
  expect_equal(m_fx@fitmeasures$df, unname(lavaan::fitMeasures(fitL)["df"]))
  expect_equal(m_fx@fitmeasures$npar, unname(lavaan::fitMeasures(fitL)["npar"]))
  expect_true(abs(m_fx@fitmeasures$logl - unname(lavaan::fitMeasures(fitL)["logl"])) < 1e-4)
  expect_true(abs(m_fx@fitmeasures$chisq - unname(lavaan::fitMeasures(fitL)["chisq"])) < 1e-4)

  # (2) lvm SEM with latent regressed on observed exogenous covariates: exact
  # estimates, SEs, chisq, df, logl.
  fitS <- lavaan::sem("dem60 =~ y1+y2+y3+y4\ndem60 ~ x1+x2+x3",
                      data = dat2, fixed.x = TRUE, meanstructure = TRUE)
  fmS <- lavaan::fitMeasures(fitS)
  peS <- lavaan::parameterEstimates(fitS)
  expect_equal(m_lvm_fx@fitmeasures$df, unname(fmS["df"]))
  expect_equal(m_lvm_fx@fitmeasures$npar, unname(fmS["npar"]))
  expect_true(abs(m_lvm_fx@fitmeasures$chisq - unname(fmS["chisq"])) < 1e-2)
  expect_true(abs(m_lvm_fx@fitmeasures$logl - unname(fmS["logl"])) < 1e-2)

  # Regression coefficients (beta) and their SEs match lavaan:
  pb <- m_lvm_fx@parameters[m_lvm_fx@parameters$matrix == "beta" &
                              m_lvm_fx@parameters$par != 0, ]
  lav_reg <- peS[peS$op == "~", ]
  # Map: Eta order is dem60(1), Lx1(2), Lx2(3), Lx3(4); beta rows are dem60<-Lx*.
  for (j in seq_len(nrow(lav_reg))){
    rhs <- lav_reg$rhs[j]               # x1/x2/x3
    latname <- c(x1 = "Lx1", x2 = "Lx2", x3 = "Lx3")[[rhs]]
    pr <- pb[pb$var2 == latname, ]
    if (nrow(pr) == 1){
      expect_true(abs(pr$est - lav_reg$est[j]) < 5e-3)
      expect_true(abs(pr$se - lav_reg$se[j]) < 1e-3)
    }
  }
}

## =========================================================================
## Regression (0.16.1): the baseline model of a fixed_x fit must be
## conditioned on x (matching lavaan). Previously baseline.df was negative
## and TLI was NaN.
## =========================================================================
if (at_home() && requireNamespace("lavaan", quietly = TRUE)){
  PD <- lavaan::PoliticalDemocracy
  m_fxb <- suppressWarnings(runmodel(
    varcov(PD, vars = c("y1","x1","x2","x3"), type = "cov",
           fixed_x = c("x1","x2","x3")), verbose = FALSE))
  l_fxb <- lavaan::sem("y1 ~ x1 + x2 + x3", PD, fixed.x = TRUE, meanstructure = TRUE)
  lfm <- lavaan::fitMeasures(l_fxb)
  expect_equal(m_fxb@fitmeasures$baseline.df, as.integer(lfm["baseline.df"]))
  expect_true(m_fxb@fitmeasures$baseline.df > 0)
  expect_true(abs(m_fxb@fitmeasures$baseline.chisq -
                    as.numeric(lfm["baseline.chisq"])) < 1e-2)
  expect_true(is.finite(m_fxb@fitmeasures$cfi))
}

## =========================================================================
## Regression (0.16.1, V2-5): an lvm fixed_x model with regressions among
## endogenous latents matches lavaan sem(fixed.x=TRUE) EXACTLY once the
## endogenous residual covariance is freed via sigma_zeta (psychonetrics,
## unlike lavaan's auto.cov.lv.y, does not free it by default).
## =========================================================================
if (at_home() && requireNamespace("lavaan", quietly = TRUE)){
  HS <- lavaan::HolzingerSwineford1939
  Lam <- matrix(0,8,4); Lam[1:3,1]<-1; Lam[4:6,2]<-1; Lam[7,3]<-1; Lam[8,4]<-1
  Bet <- matrix(0,4,4); Bet[1,3]<-Bet[1,4]<-Bet[2,3]<-Bet[2,4]<-1
  Se  <- diag(8); Se[7,7]<-0; Se[8,8]<-0
  Sz  <- matrix(0,4,4); Sz[1,1]<-Sz[2,2]<-1; Sz[2,1]<-Sz[1,2]<-1; Sz[3,3]<-Sz[4,4]<-1; Sz[4,3]<-Sz[3,4]<-1
  m <- suppressWarnings(runmodel(lvm(HS, lambda=Lam, vars=paste0("x",1:8), beta=Bet,
        latents=c("f1","f2","x7l","x8l"), sigma_epsilon=Se, sigma_zeta=Sz,
        identification="loadings", fixed_x=c("x7l","x8l")), verbose=FALSE))
  l <- lavaan::sem("f1 =~ x1+x2+x3\nf2 =~ x4+x5+x6\nf1 ~ x7+x8\nf2 ~ x7+x8",
                   HS, fixed.x=TRUE, meanstructure=TRUE)
  lfm <- lavaan::fitMeasures(l)
  expect_equal(m@fitmeasures$df, as.integer(lfm["df"]))
  expect_true(abs(m@fitmeasures$chisq - as.numeric(lfm["chisq"])) < 1e-2)
  expect_true(abs(m@fitmeasures$cfi   - as.numeric(lfm["cfi"]))   < 1e-3)
  pb <- sort(m@parameters$est[m@parameters$matrix=="beta" & !m@parameters$fixed])
  lb <- sort(lavaan::parameterEstimates(l)$est[lavaan::parameterEstimates(l)$op=="~"])
  expect_true(max(abs(pb - lb)) < 1e-3)
}
