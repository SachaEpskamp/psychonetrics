# Tests for beta_epsilon (residual temporal effects) in dlvm1, added in
# 0.16.3: the default "zero" must reproduce the white-noise residual model
# exactly; "diag" adds one residual autoregression per indicator.

suppressMessages(library(psychonetrics))

## ---- Simulate a small panel with item-level AR residuals ----
set.seed(101)
n <- 4; q <- 2; nT <- 4; N <- 500
Lam <- matrix(0, n, q); Lam[1:2, 1] <- c(1, .8); Lam[3:4, 2] <- c(1, .7)
B <- matrix(c(.3, .1, -.05, .25), q, q)
SigW <- matrix(c(.6, .15, .15, .5), q, q)
SigZb <- matrix(c(.4, .1, .1, .3), q, q)
Su <- diag(runif(n, .3, .5))
Beps_true <- diag(c(.4, .3, .25, .35))
Seb <- diag(runif(n, .15, .25))
Sig0 <- matrix(solve(diag(q^2) - kronecker(B, B), as.vector(SigW)), q, q)
Se0 <- matrix(solve(diag(n^2) - kronecker(Beps_true, Beps_true), as.vector(Su)), n, n)
bi <- MASS::mvrnorm(N, rep(0, q), SigZb); ebi <- MASS::mvrnorm(N, rep(0, n), Seb)
eta <- MASS::mvrnorm(N, rep(0, q), Sig0); eps <- MASS::mvrnorm(N, rep(0, n), Se0)
Y <- matrix(0, N, n * nT)
for (t in 1:nT){
  if (t > 1){
    eta <- eta %*% t(B) + MASS::mvrnorm(N, rep(0, q), SigW)
    eps <- eps %*% t(Beps_true) + MASS::mvrnorm(N, rep(0, n), Su)
  }
  Y[, (t - 1) * n + 1:n] <- (eta + bi) %*% t(Lam) + ebi + eps
}
colnames(Y) <- paste0("y", rep(1:n, nT), "_", rep(1:nT, each = n))
d <- as.data.frame(Y)
design <- matrix(colnames(d), n, nT); rownames(design) <- paste0("y", 1:n)
lampat <- 1 * (Lam != 0)

## ---- Default model: beta_epsilon present, fixed to zero ----
m0 <- suppressWarnings(runmodel(dlvm1(d, vars = design, lambda = lampat), verbose = FALSE))
be_rows <- m0@parameters[m0@parameters$matrix == "beta_epsilon", ]
expect_true(nrow(be_rows) == n^2)
expect_true(all(be_rows$fixed))
expect_true(all(be_rows$est == 0))

## ---- Explicit "zero" is identical to the default ----
m0b <- suppressWarnings(runmodel(dlvm1(d, vars = design, lambda = lampat,
                                       beta_epsilon = "zero"), verbose = FALSE))
expect_equal(m0@fitmeasures$logl, m0b@fitmeasures$logl)
expect_equal(m0@fitmeasures$df, m0b@fitmeasures$df)

## ---- "diag" frees one residual autoregression per indicator ----
m1 <- suppressWarnings(runmodel(dlvm1(d, vars = design, lambda = lampat,
                                      beta_epsilon = "diag"), verbose = FALSE))
expect_equal(m1@fitmeasures$npar, m0@fitmeasures$npar + n)
expect_equal(m1@fitmeasures$df, m0@fitmeasures$df - n)
expect_true(is.finite(m1@fitmeasures$logl))
# The AR model fits better (data were generated with residual AR):
expect_true(m1@fitmeasures$logl > m0@fitmeasures$logl)
dchisq <- m0@fitmeasures$chisq - m1@fitmeasures$chisq
expect_true(pchisq(dchisq, n, lower.tail = FALSE) < 1e-6)

## ---- Stationary residual covariance satisfies the Lyapunov equation ----
st <- getmatrix(m1, "sigma_epsilon_within_stationary")
iv <- getmatrix(m1, "sigma_epsilon_within")
be <- getmatrix(m1, "beta_epsilon")
expect_true(max(abs(st - (be %*% st %*% t(be) + iv))) < 1e-8)
# Off-diagonal beta_epsilon stayed fixed at zero under "diag":
expect_true(all(be[row(be) != col(be)] == 0))

## ---- Two waves: identification warning ----
expect_warning(
  dlvm1(d[, as.vector(design[, 1:2])], vars = design[, 1:2], lambda = lampat,
        beta_epsilon = "diag", baseline_saturated = FALSE),
  pattern = "at least three waves")

## ---- at_home: recovery of the generating residual autoregressions ----
if (at_home()){
  be_hat <- diag(getmatrix(m1, "beta_epsilon"))
  expect_true(max(abs(be_hat - diag(Beps_true))) < 0.15)
  # Gradient near zero at the optimum (C++ engine):
  grad <- psychonetrics:::psychonetrics_gradient_cpp(
    psychonetrics:::parVector(m1), m1)
  expect_true(mean(abs(grad)) < 0.01)
  # R and C++ engines agree:
  mR <- suppressWarnings(runmodel(usecpp(dlvm1(d, vars = design, lambda = lampat,
                                               beta_epsilon = "diag"), FALSE), verbose = FALSE))
  expect_true(abs(mR@fitmeasures$logl - m1@fitmeasures$logl) < 1e-5)
}
