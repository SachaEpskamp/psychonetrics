# Tests analogous to psychonetrics.org example 19 (graphical VAR), using
# simulated data with fixed seeds (the original example needs an OSF
# download). Also covers panelgvar() and, at_home(), the panelvar
# convergence regression from the RI-CLPM audit (4 variables, 4 waves;
# this seed combination previously diverged).
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))

## ---- Example 19 analogue: gvar pipeline on simulated VAR(1) data ----
set.seed(42)
p <- 5; TT <- 200
B <- diag(0.3, p); B[1, 2] <- 0.2; B[3, 4] <- -0.2; B[5, 1] <- 0.15
SigZ <- diag(p) * 0.8; SigZ[1, 3] <- SigZ[3, 1] <- 0.3
Y <- matrix(0, TT, p)
Y[1, ] <- MASS::mvrnorm(1, rep(0, p), diag(p))
for (t in 2:TT) Y[t, ] <- as.vector(B %*% Y[t - 1, ]) + MASS::mvrnorm(1, rep(0, p), SigZ)
tsdata <- as.data.frame(Y)
colnames(tsdata) <- paste0("V", 1:p)

model <- suppressWarnings(runmodel(gvar(tsdata, vars = paste0("V", 1:p))))
expect_inherits(model, "psychonetrics")
expect_true(abs(model@fitmeasures$logl - (-2623.75131)) < 0.01)
expect_equal(model@fitmeasures$npar, 65)
expect_equal(model@fitmeasures$df, 0)

model_pruned <- suppressWarnings(prune(model, alpha = 0.01))
expect_equal(model_pruned@fitmeasures$df, 29)

temporal <- getmatrix(model_pruned, "PDC")
contemporaneous <- getmatrix(model_pruned, "omega_zeta")
expect_equal(dim(temporal), c(5L, 5L))
expect_equal(dim(contemporaneous), c(5L, 5L))
expect_true(all(diag(contemporaneous) == 0))

cmp <- compare(saturated = model, pruned = model_pruned)
expect_true(all(is.finite(cmp$AIC)))
expect_true(all(is.finite(cmp$BIC)))

## ---- panelgvar on simulated random-intercept panel data ----
set.seed(1)
p2 <- 3; nT <- 4; N <- 300
B2 <- matrix(c(0.3, 0.1, 0,  0, 0.25, 0.1,  0.1, 0, 0.3), p2, p2, byrow = TRUE)
SigZ2 <- diag(p2) * 0.7; SigZ2[1, 2] <- SigZ2[2, 1] <- 0.2
SigB2 <- diag(p2) * 0.4
Sig02 <- matrix(solve(diag(p2^2) - kronecker(B2, B2), as.vector(SigZ2)), p2, p2)
RI <- MASS::mvrnorm(N, rep(0, p2), SigB2)
out <- matrix(NA_real_, N, p2 * nT)
for (i in 1:N) {
  eta <- MASS::mvrnorm(1, rep(0, p2), Sig02)
  out[i, 1:p2] <- eta + RI[i, ]
  for (t in 2:nT) {
    eta <- as.vector(B2 %*% eta) + MASS::mvrnorm(1, rep(0, p2), SigZ2)
    out[i, (t - 1) * p2 + 1:p2] <- eta + RI[i, ]
  }
}
colnames(out) <- paste0("V", rep(1:p2, nT), "_", rep(1:nT, each = p2))
dpanel <- as.data.frame(out)
design <- matrix(colnames(dpanel), nrow = p2, ncol = nT)
rownames(design) <- paste0("V", 1:p2)

mod_pg <- suppressWarnings(runmodel(panelgvar(dpanel, vars = design)))
expect_inherits(mod_pg, "psychonetrics")
# Since 0.16.2 panelgvar/panelvar use their own framework (not dlvm1 with
# dummy matrices); the reference values below were obtained via the dlvm1
# route and must reproduce exactly:
expect_equal(mod_pg@model, "panelvar")
expect_true(abs(mod_pg@fitmeasures$logl - (-4981.73663)) < 0.01)
expect_equal(mod_pg@fitmeasures$df, 66)
expect_equal(mod_pg@fitmeasures$npar, 24)
# The temporal network recovers the generating coefficients:
beta_hat <- getmatrix(mod_pg, "beta")
expect_equal(dim(beta_hat), c(3L, 3L))
expect_true(max(abs(beta_hat - B2)) < 0.15)
# The observed stationary means are now in 'mu' (previously 'nu'):
expect_true(all(mod_pg@parameters$matrix[mod_pg@parameters$op == "~1"] == "mu"))
# PDC is available from the implied matrices:
expect_equal(dim(getmatrix(mod_pg, "PDC")), c(3L, 3L))

## ---- at_home: exact equivalence with the old dlvm1 dummy-matrix route ----
if (at_home()) {
  I3 <- diag(p2); O3 <- matrix(0, p2, p2)
  mod_dlvm1 <- suppressWarnings(runmodel(dlvm1(
    dpanel, vars = design, lambda = I3,
    within_latent = "ggm", within_residual = "cov", sigma_epsilon_within = O3,
    between_latent = "ggm", between_residual = "cov", sigma_epsilon_between = O3
  )))
  expect_true(abs(mod_pg@fitmeasures$logl - mod_dlvm1@fitmeasures$logl) < 1e-6)
  expect_equal(mod_pg@fitmeasures$df, mod_dlvm1@fitmeasures$df)
  expect_equal(mod_pg@fitmeasures$npar, mod_dlvm1@fitmeasures$npar)
  expect_true(max(abs(getmatrix(mod_pg, "beta") - getmatrix(mod_dlvm1, "beta"))) < 1e-6)
  expect_true(max(abs(getmatrix(mod_pg, "omega_zeta_within") -
                        getmatrix(mod_dlvm1, "omega_zeta_within"))) < 1e-6)
  expect_true(max(abs(getmatrix(mod_pg, "omega_zeta_between") -
                        getmatrix(mod_dlvm1, "omega_zeta_between"))) < 1e-6)
}

## ---- at_home: panelvar convergence regression (RI-CLPM audit seeds) ----
if (at_home()) {
  # 4-variable, 4-wave panel VAR; with these seeds the model previously
  # diverged (objective stuck, huge gradient). It must now converge.
  set.seed(10)
  p3 <- 4; nT3 <- 4; N3 <- 300
  B3 <- matrix(c(0.3, 0.1, 0, 0,
                 0,   0.25, 0.1, 0,
                 0.1, 0,   0.3, -0.1,
                 0,   0.1, 0,   0.2), p3, p3, byrow = TRUE)
  Pz <- matrix(rnorm(p3 * p3), p3, p3); SigZ3 <- crossprod(Pz) / p3 + diag(p3) * 0.4
  Pb <- matrix(rnorm(p3 * p3), p3, p3); SigB3 <- crossprod(Pb) / (2 * p3) + diag(p3) * 0.3
  Sig03 <- matrix(solve(diag(p3^2) - kronecker(B3, B3), as.vector(SigZ3)), p3, p3)
  set.seed(11)
  RI3 <- MASS::mvrnorm(N3, rep(0, p3), SigB3)
  out3 <- matrix(NA_real_, N3, p3 * nT3)
  for (i in 1:N3) {
    eta <- MASS::mvrnorm(1, rep(0, p3), Sig03)
    out3[i, 1:p3] <- eta + RI3[i, ]
    for (t in 2:nT3) {
      eta <- as.vector(B3 %*% eta) + MASS::mvrnorm(1, rep(0, p3), SigZ3)
      out3[i, (t - 1) * p3 + 1:p3] <- eta + RI3[i, ]
    }
  }
  colnames(out3) <- paste0("V", rep(1:p3, nT3), "_", rep(1:nT3, each = p3))
  d3 <- as.data.frame(out3)
  design3 <- matrix(colnames(d3), nrow = p3, ncol = nT3)
  rownames(design3) <- paste0("V", 1:p3)

  mod_pv <- suppressWarnings(runmodel(panelvar(d3, vars = design3)))
  ll <- mod_pv@fitmeasures$logl
  expect_true(is.finite(ll))
  expect_true(ll < 0)
  expect_true(abs(ll - (-7106.99383)) < 1)
  grad <- psychonetrics:::psychonetrics_gradient_cpp(
    psychonetrics:::parVector(mod_pv), mod_pv
  )
  expect_true(mean(abs(grad)) < 0.01)
}
