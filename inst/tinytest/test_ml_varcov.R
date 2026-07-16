# Tests for the multi-level variance-covariance family (ml_varcov / ml_ggm /
# ml_corr, experimental). ml_varcov is the multi-level analogue of the varcov
# family, estimated with the ml_lvm two-level sufficient-statistics ML
# machinery. Fast deterministic checks run always (CRAN); the heavier loops over
# all type combinations and the ml_lvm cross-check run under at_home().

suppressMessages({
  library(psychonetrics)
})

# ---------------------------------------------------------------------------
# Simulate two-level data with a known within/between structure (complete):
# ---------------------------------------------------------------------------
sim_ml <- function(seed, p, J, nper, SigW, SigB, mu = rep(0, p)){
  set.seed(seed)
  cl <- rep(seq_len(J), each = nper)
  b  <- MASS::mvrnorm(J, rep(0, p), SigB)
  Y  <- b[cl, , drop = FALSE] + MASS::mvrnorm(J * nper, rep(0, p), SigW)
  Y  <- sweep(Y, 2, mu, "+")
  d  <- as.data.frame(Y); names(d) <- paste0("V", seq_len(p)); d$cluster <- cl
  d
}

p <- 4L; J <- 200L; nper <- 6L
Ow <- matrix(0, p, p); Ow[1,2] <- Ow[2,1] <- .3; Ow[2,3] <- Ow[3,2] <- .25; Ow[3,4] <- Ow[4,3] <- -.2
SigW0 <- cov2cor(solve(diag(p) - Ow))
SigB0 <- diag(c(.5, .6, .4, .5)); SigB0[1,2] <- SigB0[2,1] <- .2
dat <- sim_ml(1, p, J, nper, SigW0, SigB0, mu = c(1, 2, 3, 4))

# Column-major vech (lower triangle incl. diagonal):
vechf <- function(M) M[lower.tri(M, diag = TRUE)]

# ---------------------------------------------------------------------------
# (1) Construction: model/submodel/type slots and native matrix names
# ---------------------------------------------------------------------------
mod <- ml_ggm(dat, vars = paste0("V", 1:p), clusters = "cluster")
expect_equal(mod@model, "ml_varcov")
expect_equal(mod@submodel, "ml_ggm")
expect_equal(mod@types$within, "ggm")
expect_equal(mod@types$between, "ggm")
# Native names, no "zeta":
mnames <- unique(mod@parameters$matrix)
expect_true(all(c("mu", "omega_within", "delta_within", "omega_between", "delta_between") %in% mnames))
expect_false(any(grepl("zeta|epsilon|lambda|beta", mnames)))

# ml_corr wrapper sets type "cor":
modc <- ml_corr(dat, vars = paste0("V", 1:p), clusters = "cluster")
expect_equal(modc@submodel, "ml_corr")
expect_equal(modc@types$within, "cor")

# ---------------------------------------------------------------------------
# (2) Saturated ml_ggm fits with 0 df and recovers the sufficient statistics
# ---------------------------------------------------------------------------
mod <- runmodel(mod, verbose = FALSE)
expect_true(mod@computed)
expect_equal(fit(mod)$Value[fit(mod)$Measure == "df"], 0)

# ---------------------------------------------------------------------------
# (3) Analytic gradient == numeric gradient (R path), at start and a random point
# ---------------------------------------------------------------------------
grad_check <- function(m, tol = 1e-5){
  m <- usecpp(m, FALSE)
  xs <- psychonetrics:::parVector(m)
  gA <- psychonetrics:::psychonetrics_gradient(xs, m)
  gN <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, m), xs)
  ok1 <- max(abs(gA - gN)) < tol
  set.seed(7); xr <- xs + rnorm(length(xs), 0, 0.03)
  gAr <- psychonetrics:::psychonetrics_gradient(xr, m)
  gNr <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, m), xr)
  ok2 <- max(abs(gAr - gNr)) < tol
  ok1 && ok2
}
expect_true(grad_check(mod))

# A restricted (non-saturated) model gives finite, sane fit and correct gradient:
Wpat <- matrix(1, p, p); Wpat[1,4] <- Wpat[4,1] <- 0; Wpat[1,3] <- Wpat[3,1] <- 0; diag(Wpat) <- 0
Bpat <- diag(p); Bpat[1,2] <- Bpat[2,1] <- 1
modR <- ml_varcov(dat, vars = paste0("V", 1:p), clusters = "cluster",
                  within = "ggm", omega_within = Wpat,
                  between = "ggm", omega_between = Bpat)
modR <- runmodel(modR, verbose = FALSE)
expect_true(fit(modR)$Value[fit(modR)$Measure == "df"] > 0)
expect_true(grad_check(modR))

# ---------------------------------------------------------------------------
# (4) R vs C++ agreement (only when a C++ build is available)
# ---------------------------------------------------------------------------
mod_cpp <- usecpp(mod, TRUE)
if (isTRUE(mod_cpp@cpp)){
  xs <- psychonetrics:::parVector(mod_cpp)
  gR <- psychonetrics:::psychonetrics_gradient(xs, usecpp(mod, FALSE))
  gC <- psychonetrics:::psychonetrics_gradient_cpp(xs, mod_cpp)
  expect_true(max(abs(gR - gC)) < 1e-8)
}

# ---------------------------------------------------------------------------
# (5) FIML estimator: construction, gradient, and exact equivalence with
#     ml_lvm-FIML at identical parameter vectors
# ---------------------------------------------------------------------------
modF <- ml_ggm(dat, vars = paste0("V", 1:p), clusters = "cluster",
               estimator = "FIML", baseline_saturated = FALSE)
expect_equal(modF@estimator, "FIML")
expect_equal(modF@distribution, "Gaussian")
expect_true(grad_check(modF))

# The decisive equivalence check: ml_varcov-FIML and ml_lvm-FIML (with lambda
# missing, i.e. the identity) are the SAME model with aligned parameter
# vectors, so their fit functions and gradients must agree to machine
# precision at any theta:
modL <- usecpp(suppressMessages(ml_lvm(dat, vars = paste0("V", 1:p), clusters = "cluster",
                     within_latent = "ggm", between_latent = "ggm",
                     estimator = "FIML", baseline_saturated = FALSE)), FALSE)
modF_R <- usecpp(modF, FALSE)
xs <- psychonetrics:::parVector(modF_R)
expect_equal(length(xs), length(psychonetrics:::parVector(modL)))
set.seed(3); xr <- xs + rnorm(length(xs), 0, 0.04)
fF <- psychonetrics:::psychonetrics_fitfunction(xr, modF_R)
fL <- psychonetrics:::psychonetrics_fitfunction(xr, modL)
expect_true(abs(fF - fL) < 1e-10)
gF <- psychonetrics:::psychonetrics_gradient(xr, modF_R)
gL <- psychonetrics:::psychonetrics_gradient(xr, modL)
expect_true(max(abs(gF - gL)) < 1e-10)

# ---------------------------------------------------------------------------
# (6) Row order and cluster labels: CLUSTERID regression tests. Shuffling the
#     rows must not change the model (this crashed before the ave() fix), and
#     cluster labels reused across groups must denote different clusters:
# ---------------------------------------------------------------------------
set.seed(4); dat_shuf <- dat[sample(nrow(dat)), ]
modS <- ml_ggm(dat_shuf, vars = paste0("V", 1:p), clusters = "cluster",
               estimator = "FIML", baseline_saturated = FALSE)
fS <- psychonetrics:::psychonetrics_fitfunction(xr, usecpp(modS, FALSE))
expect_true(abs(fS - fF) < 1e-10)

dat2g <- rbind(cbind(dat, g = "A"), cbind(dat, g = "B"))  # reused cluster labels
mod2g <- ml_ggm(dat2g, vars = paste0("V", 1:p), clusters = "cluster",
                groupvar = "g", estimator = "FIML", baseline_saturated = FALSE)
expect_equal(nrow(mod2g@sample@groups), 2L)
expect_equal(mod2g@sample@groups$nobs, c(J, J))  # J clusters per group, not shared

# ---------------------------------------------------------------------------
# Heavier checks (all type combinations; ml_lvm cross-validation; FIML
# missing data; ML/FIML agreement):
# ---------------------------------------------------------------------------
if (at_home()){
  types <- c("cov", "chol", "prec", "ggm", "cor")
  for (w in types) for (b_ in types){
    m <- ml_varcov(dat, vars = paste0("V", 1:p), clusters = "cluster", estimator = "ML",
                   within = w, between = b_, baseline_saturated = FALSE)
    expect_true(grad_check(m), info = paste0("ML gradient within=", w, " between=", b_))
    mf <- ml_varcov(dat, vars = paste0("V", 1:p), clusters = "cluster", estimator = "FIML",
                    within = w, between = b_, baseline_saturated = FALSE)
    expect_true(grad_check(mf), info = paste0("FIML gradient within=", w, " between=", b_))
  }

  # ml_ggm equals ml_lvm with lambda missing (same model, different framework):
  m_vc <- runmodel(ml_ggm(dat, vars = paste0("V", 1:p), clusters = "cluster", estimator = "ML"), verbose = FALSE)
  m_lv <- runmodel(suppressMessages(ml_lvm(dat, vars = paste0("V", 1:p), clusters = "cluster",
                          within_latent = "ggm", between_latent = "ggm", estimator = "ML")), verbose = FALSE)
  expect_equal(max(abs(getmatrix(m_vc, "sigma_within")  - getmatrix(m_lv, "sigma_within"))),  0, tolerance = 1e-3)
  expect_equal(max(abs(getmatrix(m_vc, "sigma_between") - getmatrix(m_lv, "sigma_between"))), 0, tolerance = 1e-3)

  # ML and FIML land on the same optimum on complete data (same likelihood):
  m_fi <- runmodel(ml_ggm(dat, vars = paste0("V", 1:p), clusters = "cluster", estimator = "FIML"), verbose = FALSE)
  expect_equal(m_vc@fitmeasures$logl, m_fi@fitmeasures$logl, tolerance = 1e-4)
  expect_equal(max(abs(psychonetrics:::parVector(m_vc) - psychonetrics:::parVector(m_fi))), 0, tolerance = 1e-3)

  # FIML with missing data: gradient correct and equal to ml_lvm-FIML:
  datm <- dat
  set.seed(5); datm[sample(nrow(datm), round(0.07 * nrow(datm))), sample(p, 1)] <- NA
  mm <- suppressWarnings(ml_ggm(datm, vars = paste0("V", 1:p), clusters = "cluster",
                                estimator = "FIML", baseline_saturated = FALSE))
  expect_true(grad_check(mm))
  mmL <- usecpp(suppressWarnings(suppressMessages(ml_lvm(datm, vars = paste0("V", 1:p), clusters = "cluster",
              within_latent = "ggm", between_latent = "ggm",
              estimator = "FIML", baseline_saturated = FALSE))), FALSE)
  xm <- psychonetrics:::parVector(mm)
  set.seed(6); xmr <- xm + rnorm(length(xm), 0, 0.04)
  expect_true(abs(psychonetrics:::psychonetrics_fitfunction(xmr, usecpp(mm, FALSE)) -
                  psychonetrics:::psychonetrics_fitfunction(xmr, mmL)) < 1e-10)

  # prune() and stepup() run and return a psychonetrics object:
  expect_inherits(prune(m_vc, alpha = 0.01, verbose = FALSE), "psychonetrics")
}
