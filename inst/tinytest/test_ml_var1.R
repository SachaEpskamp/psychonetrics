# Tests for the multi-level lag-1 VAR / graphical VAR family (ml_var1 /
# ml_gvar1, experimental, added in 0.16.6). This is the var1 lag-embedding
# pseudo-likelihood generalized with a between-person random-intercept level,
# estimated by the ml_lvm two-level summary-statistics ML machinery.
#
# Fast deterministic checks run always (CRAN); the heavier loops over all
# type combinations and the recovery smoke test run under at_home().

suppressMessages(library(psychonetrics))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Column-major lower-triangle-including-diagonal vech:
vechf <- function(M) M[lower.tri(M, diag = TRUE)]

# Independent clean-room two-level log-likelihood from sufficient statistics
# (McDonald-Goldstein / Muthen), evaluated at (mu, SW, SB) with cluster = id.
# This is the reference implementation for the fit function (design doc 6.5).
cleanroom_logl_2l <- function(Y, cluster, mu, SW, SB){
  Y <- as.matrix(Y); N <- nrow(Y); q <- ncol(Y)
  cl <- as.integer(factor(cluster)); J <- max(cl)
  ns <- tabulate(cl)
  Ybar <- rowsum(Y, cl) / ns
  S.PW <- crossprod(Y - Ybar[cl, , drop = FALSE]) / (N - J)
  iSW <- solve(SW)
  m2ll <- (N - J) * (as.numeric(determinant(SW, TRUE)$modulus) + sum(iSW * S.PW))
  for (s in sort(unique(ns))){
    idx <- which(ns == s); m <- length(idx)
    Yb <- Ybar[idx, , drop = FALSE]
    mb <- colMeans(Yb)
    Cb <- if (m > 1) crossprod(sweep(Yb, 2, mb)) / m else matrix(0, q, q)
    Sj <- SW + s * SB
    iSj <- solve(Sj)
    A <- Cb + tcrossprod(mb - mu)
    m2ll <- m2ll + m * (as.numeric(determinant(Sj, TRUE)$modulus) + s * sum(iSj * A))
  }
  -0.5 * (m2ll + N * q * log(2 * pi))
}

# Simulate a two-level VAR(1) data set with random intercepts (long format):
sim_mlvar <- function(seed, p, J, TT, beta, SigZ, SigB, mu = rep(0, p)){
  set.seed(seed)
  Sig0 <- matrix(solve(diag(p^2) - kronecker(beta, beta), as.vector(SigZ)), p, p)
  long <- do.call(rbind, lapply(seq_len(J), function(i){
    bi  <- MASS::mvrnorm(1, rep(0, p), SigB)
    eta <- MASS::mvrnorm(1, rep(0, p), Sig0)
    Y <- matrix(NA_real_, TT, p)
    Y[1, ] <- mu + eta + bi
    for (t in 2:TT){
      eta <- as.vector(beta %*% eta) + MASS::mvrnorm(1, rep(0, p), SigZ)
      Y[t, ] <- mu + eta + bi
    }
    data.frame(id = i, beep = seq_len(TT), Y)
  }))
  names(long)[2 + seq_len(p)] <- paste0("V", seq_len(p))
  long
}

# A small fixed data set reused across tests:
p <- 3L; J <- 40L; TT <- 25L
beta0 <- diag(0.3, p); beta0[1, 2] <- 0.2; beta0[3, 1] <- -0.15
SigZ0 <- diag(p) * 0.8; SigZ0[1, 2] <- SigZ0[2, 1] <- 0.2
SigB0 <- diag(p) * 0.4; SigB0[2, 3] <- SigB0[3, 2] <- 0.1
dat <- sim_mlvar(1, p, J, TT, beta0, SigZ0, SigB0)

# ---------------------------------------------------------------------------
# (1) tsData(includeID) backward compatibility
# ---------------------------------------------------------------------------
ts_default <- psychonetrics:::tsData(dat, vars = paste0("V", 1:p),
                                     idvar = "id", beepvar = "beep")
ts_F <- psychonetrics:::tsData(dat, vars = paste0("V", 1:p),
                               idvar = "id", beepvar = "beep", includeID = FALSE)
ts_T <- psychonetrics:::tsData(dat, vars = paste0("V", 1:p),
                               idvar = "id", beepvar = "beep", includeID = TRUE)

# includeID = FALSE reproduces the historical output exactly:
expect_identical(ts_default, ts_F)
# includeID = TRUE adds exactly the id column, leaving the rest unchanged:
expect_true("id" %in% colnames(ts_T))
expect_false("id" %in% colnames(ts_F))
expect_equal(ncol(ts_T), ncol(ts_F) + 1L)
expect_equal(ts_T[, colnames(ts_F), drop = FALSE], ts_F,
             check.attributes = FALSE)

# ---------------------------------------------------------------------------
# Build a base model (construction only) to access the structure maps and
# extra matrices:
# ---------------------------------------------------------------------------
mod0 <- suppressMessages(suppressWarnings(
  ml_gvar1(dat, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep")))
expect_inherits(mod0, "psychonetrics")
expect_equal(mod0@model, "ml_var1")
expect_equal(mod0@submodel, "ml_gvar1")

# ---------------------------------------------------------------------------
# (2) Structure-map identities P_within / P_between
# ---------------------------------------------------------------------------
Pw <- as.matrix(mod0@extramatrices$P_within)
Pb <- as.matrix(mod0@extramatrices$P_between)

set.seed(123)
Araw <- matrix(rnorm(p^2), p, p); S0 <- crossprod(Araw) + diag(p)   # symmetric PD
S1 <- matrix(rnorm(p^2), p, p)                                       # non-symmetric
SW <- rbind(cbind(S0, t(S1)), cbind(S1, S0))
# vech(Sigma_W) == P_within %*% [ vech(Sigma0) ; vec(Sigma1) ]
expect_true(max(abs(as.numeric(Pw %*% c(vechf(S0), as.vector(S1))) - vechf(SW))) < 1e-10)

SB <- crossprod(matrix(rnorm(p^2), p, p)) + diag(p)
SBf <- kronecker(matrix(1, 2, 2), SB)
# vech(1_2x2 (x) Sigma_B) == P_between %*% vech(Sigma_B)
expect_true(max(abs(as.numeric(Pb %*% vechf(SB)) - vechf(SBf))) < 1e-10)

# ---------------------------------------------------------------------------
# (6) nobs equals number of subjects; df = nGroup*(2p + 2k) - npar
# ---------------------------------------------------------------------------
expect_equal(mod0@sample@groups$nobs[1], J)               # subjects, not rows
k <- p * (2 * p + 1)
npar0 <- max(mod0@parameters$par)
expect_equal(mod0@sample@nobs, 1L * (2 * p + 2 * k))

# ---------------------------------------------------------------------------
# (7) Multi-group: equal = "beta" reduces the free-parameter count by p^2
# ---------------------------------------------------------------------------
dat2 <- rbind(transform(dat, g = "A"),
              transform(sim_mlvar(2, p, J, TT, beta0, SigZ0, SigB0), g = "B",
                        id = id + 1000L))
mod_2g <- suppressMessages(suppressWarnings(
  ml_gvar1(dat2, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep",
           groupvar = "g")))
mod_2g_eq <- suppressMessages(suppressWarnings(
  ml_gvar1(dat2, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep",
           groupvar = "g", equal = "beta")))
np_free <- max(mod_2g@parameters$par)
np_eq   <- max(mod_2g_eq@parameters$par)
expect_true(np_eq < np_free)
expect_equal(np_free - np_eq, p^2)   # one shared beta instead of two

# ---------------------------------------------------------------------------
# (4) Analytic gradient equals numDeriv at start and at a perturbed point,
#     for a couple of type combinations.
# ---------------------------------------------------------------------------
grad_check <- function(mod, tol = 1e-5){
  xs <- psychonetrics:::parVector(mod)
  gA <- psychonetrics:::psychonetrics_gradient(xs, mod)
  gN <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, mod), xs)
  ok1 <- max(abs(gA - gN)) < tol
  set.seed(7)
  xr <- xs + rnorm(length(xs), 0, 0.03)
  gAr <- psychonetrics:::psychonetrics_gradient(xr, mod)
  gNr <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, mod), xr)
  ok2 <- max(abs(gAr - gNr)) < tol
  ok1 && ok2
}
mod_cov <- suppressMessages(suppressWarnings(
  ml_var1(dat, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep")))
expect_true(grad_check(mod_cov))   # cov / cov
expect_true(grad_check(mod0))      # ggm / ggm

# ---------------------------------------------------------------------------
# (5) Fit function value equals the independent hand-rolled reference.
#     Compare at the (converged) solution using the implied (mu, SW, SB).
# ---------------------------------------------------------------------------
mod_run <- suppressMessages(suppressWarnings(runmodel(mod_cov)))
expect_inherits(mod_run, "psychonetrics")
expect_true(is.finite(mod_run@fitmeasures$logl))
expect_true(is.finite(mod_run@fitmeasures$aic.ll))
expect_true(is.finite(mod_run@fitmeasures$bic))

# Rebuild the exact augmented + complete-pairs data used internally:
aug <- psychonetrics:::tsData(dat, vars = paste0("V", 1:p),
                              idvar = "id", beepvar = "beep", includeID = TRUE)
vars2p <- setdiff(colnames(aug), "id")
aug <- aug[rowSums(is.na(aug[, vars2p, drop = FALSE])) == 0, , drop = FALSE]

# Implied (mu, Sigma_W, Sigma_B) at the fitted point:
beta_hat <- getmatrix(mod_run, "beta")
SigZ_hat <- getmatrix(mod_run, "sigma_zeta_within")
SigB_hat <- getmatrix(mod_run, "sigma_zeta_between")
# getmatrix(..., "mu") already returns the implied 2p augmented mean:
muz_hat  <- as.numeric(getmatrix(mod_run, "mu"))
Sig0_hat <- matrix(solve(diag(p^2) - kronecker(beta_hat, beta_hat),
                         as.vector(SigZ_hat)), p, p)
Sig1_hat <- beta_hat %*% Sig0_hat
SW_hat <- rbind(cbind(Sig0_hat, t(Sig1_hat)), cbind(Sig1_hat, Sig0_hat))
SB_hat <- kronecker(matrix(1, 2, 2), SigB_hat)

ll_ref <- cleanroom_logl_2l(aug[, vars2p], aug[["id"]], muz_hat, SW_hat, SB_hat)
expect_true(abs(ll_ref - mod_run@fitmeasures$logl) < 1e-4)

# ---------------------------------------------------------------------------
# (8) Loose parameter recovery on the fitted ggm / ggm model
# ---------------------------------------------------------------------------
mod0_run <- suppressMessages(suppressWarnings(runmodel(mod0)))
beta_r <- getmatrix(mod0_run, "beta")
oW_r   <- getmatrix(mod0_run, "omega_zeta_within")
oB_r   <- getmatrix(mod0_run, "omega_zeta_between")
expect_equal(dim(beta_r), c(p, p))
expect_equal(dim(oW_r), c(p, p))
expect_equal(dim(oB_r), c(p, p))
# Temporal effects recovered with the right sign/magnitude (loose):
expect_true(abs(beta_r[1, 2] - 0.2) < 0.15)
expect_true(beta_r[3, 1] < 0)
expect_true(all(diag(oW_r) == 0) && all(diag(oB_r) == 0))

# The experimental-message mechanism exists and the model still builds:
expect_true(is.function(psychonetrics:::experimentalWarning))

# ---------------------------------------------------------------------------
# (3) Build (and, once, run) across all 5 within types x {cov, ggm} between.
# ---------------------------------------------------------------------------
if (at_home()){
  wtypes <- c("cov", "chol", "prec", "ggm", "cor")
  btypes <- c("cov", "ggm")
  for (wt in wtypes) for (bt in btypes){
    m <- suppressMessages(suppressWarnings(
      ml_var1(dat, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep",
              within_latent = wt, between_latent = bt)))
    expect_inherits(m, "psychonetrics")
    expect_true(max(m@parameters$par) > 0)
  }

  # A full pipeline for the chol/cov combination must converge:
  m_run <- suppressMessages(suppressWarnings(runmodel(
    ml_var1(dat, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep",
            within_latent = "chol", between_latent = "cov"))))
  expect_true(is.finite(m_run@fitmeasures$logl))
  gr <- psychonetrics:::psychonetrics_gradient(
    psychonetrics:::parVector(m_run), m_run)
  expect_true(mean(abs(gr)) < 0.01)
}
