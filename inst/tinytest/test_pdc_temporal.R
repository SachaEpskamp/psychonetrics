# Tests for the PDC temporal parameterization (0.16.4): temporal = "PDC"
# models the partial directed correlations (from = row, to = column)
# directly. A saturated PDC model is fit-equivalent to the raw model, and
# the fitted PDC parameters equal the derived PDC of the raw fit.

suppressMessages(library(psychonetrics))

## ---- var1: PDC vs raw equivalence ----
set.seed(42)
p <- 3; TT <- 300
B <- diag(0.3, p); B[1, 2] <- 0.2; B[3, 1] <- 0.15
SigZ <- diag(p) * 0.8; SigZ[1, 3] <- SigZ[3, 1] <- 0.3
Y <- matrix(0, TT, p)
for (t in 2:TT) Y[t, ] <- as.vector(B %*% Y[t - 1, ]) + MASS::mvrnorm(1, rep(0, p), SigZ)
tsdat <- as.data.frame(Y); colnames(tsdat) <- paste0("V", 1:p)

m_raw <- suppressWarnings(runmodel(var1(tsdat, vars = paste0("V", 1:p)), verbose = FALSE))
m_pdc <- suppressWarnings(runmodel(var1(tsdat, vars = paste0("V", 1:p), temporal = "PDC"), verbose = FALSE))

expect_equal(m_raw@fitmeasures$df, m_pdc@fitmeasures$df)
expect_equal(m_raw@fitmeasures$npar, m_pdc@fitmeasures$npar)
expect_true(abs(m_raw@fitmeasures$logl - m_pdc@fitmeasures$logl) < 1e-4)
# Fitted PDC parameters equal the derived PDC of the raw fit:
expect_true(max(abs(getmatrix(m_pdc, "PDC") - getmatrix(m_raw, "PDC"))) < 1e-3)
# beta is recoverable from the PDC model:
expect_true(max(abs(getmatrix(m_pdc, "beta") - getmatrix(m_raw, "beta"))) < 1e-3)

# PDC parameters are bounded in (-1, 1):
pp <- m_pdc@parameters[m_pdc@parameters$matrix == "PDC" & !m_pdc@parameters$fixed, ]
expect_true(all(pp$minimum == -1) && all(pp$maximum == 1))
expect_true(all(abs(pp$est) < 1))

# Parameter table orientation: from = row (var1) -> to = column (var2), op "->":
expect_true(all(pp$op == "->"))

## ---- CIplot arrow direction ----
if (requireNamespace("ggplot2", quietly = TRUE)){
  plt <- suppressWarnings(CIplot(m_pdc, matrices = "PDC", print = FALSE))
  edges_plot <- sort(as.character(unique(plt$data$edge)))
  pp_all <- m_pdc@parameters[m_pdc@parameters$matrix == "PDC", ]
  expect_equal(edges_plot, sort(unique(paste0(pp_all$var1, " -> ", pp_all$var2))))
}

## ---- transmod does not transform temporal types ----
expect_error(transmod(m_pdc, temporal = "raw"), pattern = "not supported")

## ---- at_home: panelvar and dlvm1 (latent + residual) PDC ----
if (at_home()){
  set.seed(1)
  p2 <- 3; nT <- 4; N <- 300
  B2 <- matrix(c(0.3, 0.1, 0, 0, 0.25, 0.1, 0.1, 0, 0.3), p2, p2, byrow = TRUE)
  SigZ2 <- diag(p2) * 0.7; SigZ2[1, 2] <- SigZ2[2, 1] <- 0.2
  SigB2 <- diag(p2) * 0.4
  Sig02 <- matrix(solve(diag(p2^2) - kronecker(B2, B2), as.vector(SigZ2)), p2, p2)
  RI <- MASS::mvrnorm(N, rep(0, p2), SigB2)
  out <- matrix(NA_real_, N, p2 * nT)
  for (i in 1:N){
    eta <- MASS::mvrnorm(1, rep(0, p2), Sig02)
    out[i, 1:p2] <- eta + RI[i, ]
    for (t in 2:nT){
      eta <- as.vector(B2 %*% eta) + MASS::mvrnorm(1, rep(0, p2), SigZ2)
      out[i, (t - 1) * p2 + 1:p2] <- eta + RI[i, ]
    }
  }
  colnames(out) <- paste0("V", rep(1:p2, nT), "_", rep(1:nT, each = p2))
  dpanel <- as.data.frame(out)
  design <- matrix(colnames(dpanel), p2, nT); rownames(design) <- paste0("V", 1:p2)

  mp_raw <- suppressWarnings(runmodel(panelgvar(dpanel, vars = design), verbose = FALSE))
  mp_pdc <- suppressWarnings(runmodel(panelgvar(dpanel, vars = design, temporal = "PDC"), verbose = FALSE))
  expect_true(abs(mp_raw@fitmeasures$logl - mp_pdc@fitmeasures$logl) < 1e-4)
  expect_true(max(abs(getmatrix(mp_pdc, "PDC") - getmatrix(mp_raw, "PDC"))) < 1e-3)
  # prune() defaults select the PDC matrix for PDC models:
  mp_pruned <- suppressWarnings(prune(mp_pdc, alpha = 0.01, verbose = FALSE))
  expect_true(sum(mp_pruned@parameters$fixed[mp_pruned@parameters$matrix == "PDC"]) >=
                sum(mp_pdc@parameters$fixed[mp_pdc@parameters$matrix == "PDC"]))

  # dlvm1 residual PDC equals raw beta_epsilon = "diag" via the derived
  # PDC_epsilon matrix:
  n4 <- 4; q4 <- 2
  Lam4 <- matrix(0, n4, q4); Lam4[1:2, 1] <- c(1, .8); Lam4[3:4, 2] <- c(1, .7)
  Beps <- diag(c(.4, .3, .25, .35)); Su <- diag(runif(n4, .3, .5))
  Seb <- diag(runif(n4, .15, .25)); SigW <- matrix(c(.6, .15, .15, .5), q4, q4)
  SigZb <- matrix(c(.4, .1, .1, .3), q4, q4)
  B4 <- matrix(c(.3, .1, -.05, .25), q4, q4)
  Sig04 <- matrix(solve(diag(q4^2) - kronecker(B4, B4), as.vector(SigW)), q4, q4)
  Se04 <- matrix(solve(diag(n4^2) - kronecker(Beps, Beps), as.vector(Su)), n4, n4)
  set.seed(7)
  bi <- MASS::mvrnorm(N, rep(0, q4), SigZb); ebi <- MASS::mvrnorm(N, rep(0, n4), Seb)
  eta <- MASS::mvrnorm(N, rep(0, q4), Sig04); eps <- MASS::mvrnorm(N, rep(0, n4), Se04)
  Y4 <- matrix(0, N, n4 * nT)
  for (t in 1:nT){
    if (t > 1){
      eta <- eta %*% t(B4) + MASS::mvrnorm(N, rep(0, q4), SigW)
      eps <- eps %*% t(Beps) + MASS::mvrnorm(N, rep(0, n4), Su)
    }
    Y4[, (t - 1) * n4 + 1:n4] <- (eta + bi) %*% t(Lam4) + ebi + eps
  }
  colnames(Y4) <- paste0("y", rep(1:n4, nT), "_", rep(1:nT, each = n4))
  d4 <- as.data.frame(Y4)
  des4 <- matrix(colnames(d4), n4, nT); rownames(des4) <- paste0("y", 1:n4)
  lampat4 <- 1 * (Lam4 != 0)

  me_raw <- suppressWarnings(runmodel(dlvm1(d4, vars = des4, lambda = lampat4,
                                            beta_epsilon = "diag"), verbose = FALSE))
  me_pdc <- suppressWarnings(runmodel(dlvm1(d4, vars = des4, lambda = lampat4,
                                            temporal_residual = "PDC",
                                            PDC_epsilon = "diag"), verbose = FALSE))
  expect_true(abs(me_raw@fitmeasures$logl - me_pdc@fitmeasures$logl) < 1e-3)
  expect_true(max(abs(getmatrix(me_pdc, "PDC_epsilon") - getmatrix(me_raw, "PDC_epsilon"))) < 5e-3)
}
