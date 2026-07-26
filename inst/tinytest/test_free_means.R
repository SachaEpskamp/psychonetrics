# Tests for the wave-specific ("free") mean structure added in 0.17.1:
# stationary_means = FALSE in the panelvar/panelgvar and
# dlvm1 frameworks. Covers the degrees-of-freedom bookkeeping, the saturated
# mean structure, the baseline model mirroring the target's mean structure,
# the equivalence with per-wave mean-centering, the cross-framework
# equivalence between panelgvar and dlvm1 with lambda = I, and the agreement
# between the R and the C++ estimation paths.

suppressMessages(library(psychonetrics))

# Shared helpers (portable multivariate normal draws; see _helpers.R):
source(if (file.exists("_helpers.R")) "_helpers.R" else
         system.file("tinytest", "_helpers.R", package = "psychonetrics"))

## ---- Simulate a stationary panel VAR with a random intercept and
##      wave-specific level shifts ----
sim_panel <- function(seed = 1, m = 3, nT = 3, N = 300, shifts = rep(0, nT)){
  set.seed(seed)
  B <- diag(0.3, m); B[1, 2] <- 0.15; B[2, 3] <- 0.1
  SigZ <- diag(m) * 0.6; SigZ[1, 2] <- SigZ[2, 1] <- 0.2
  SigB <- diag(m) * 0.4; SigB[1, 2] <- SigB[2, 1] <- 0.15
  Sig0 <- matrix(solve(diag(m^2) - kronecker(B, B), as.vector(SigZ)), m, m)
  RI  <- rmvnorm(N, rep(0, m), SigB)
  eta <- rmvnorm(N, rep(0, m), Sig0)
  Y <- matrix(NA_real_, N, m * nT)
  for (t in seq_len(nT)){
    if (t > 1) eta <- eta %*% t(B) + rmvnorm(N, rep(0, m), SigZ)
    Y[, (t - 1) * m + seq_len(m)] <- eta + RI + shifts[t]
  }
  colnames(Y) <- paste0("V", rep(seq_len(m), nT), "_", rep(seq_len(nT), each = m))
  d <- as.data.frame(Y)
  design <- matrix(colnames(d), m, nT)
  rownames(design) <- paste0("V", seq_len(m))
  colnames(design) <- paste0("T", seq_len(nT))
  list(data = d, design = design, m = m, nT = nT)
}

s <- sim_panel(shifts = c(0, 0.5, 1))
m <- s$m; nT <- s$nT

## ---- Default is unchanged: stationary means ----
expect_true(isTRUE(eval(formals(panelvar)$stationary_means)))
expect_true(isTRUE(eval(formals(dlvm1)$stationary_means)))

mod_stat <- suppressWarnings(runmodel(panelgvar(s$data, vars = s$design), verbose = FALSE))
mod_free <- suppressWarnings(runmodel(panelgvar(s$data, vars = s$design, stationary_means = FALSE),
                                      verbose = FALSE))

## ---- The stationary mean structure keeps a single mu column ----
expect_equal(unique(mod_stat@matrices$ncol[mod_stat@matrices$name == "mu"]), 1L)
expect_equal(unique(mod_free@matrices$ncol[mod_free@matrices$name == "mu"]), as.integer(nT))

## ---- df bookkeeping: m(T-1) extra parameters ----
expect_equal(mod_free@fitmeasures$npar, mod_stat@fitmeasures$npar + m * (nT - 1))
expect_equal(mod_stat@fitmeasures$df - mod_free@fitmeasures$df, m * (nT - 1))

## ---- 'means' keeps its meaning (sample means): a covs + means fit with
##      stationary_means = FALSE reproduces the raw-data fit ----
cv <- cov(s$data) * (nrow(s$data) - 1) / nrow(s$data)
mod_covs <- suppressWarnings(runmodel(panelgvar(vars = s$design, covs = cv,
                                                means = colMeans(s$data), nobs = nrow(s$data),
                                                covtype = "ML", stationary_means = FALSE),
                                      verbose = FALSE))
expect_true(abs(mod_covs@fitmeasures$chisq - mod_free@fitmeasures$chisq) < 1e-6)
expect_true(max(abs(mod_covs@parameters$est - mod_free@parameters$est)) < 1e-6)

## ---- A character value in 'means', and a non-logical stationary_means, error ----
expect_error(panelgvar(s$data, vars = s$design, means = "free"))
expect_error(dlvm1(s$data, vars = s$design, lambda = diag(m), means = "free"))
expect_error(panelgvar(s$data, vars = s$design, stationary_means = "free"))
expect_error(panelgvar(s$data, vars = s$design, stationary_means = c(TRUE, FALSE)))

## ---- The mean structure is exactly saturated: implied means = observed
##      per-wave means, and there is one mu parameter per variable-wave ----
mu_par <- mod_free@parameters[mod_free@parameters$matrix == "mu", ]
expect_equal(nrow(mu_par), m * nT)
expect_true(all(!mu_par$fixed))
expect_true(max(abs(as.vector(getmatrix(mod_free, "mu")) - colMeans(s$data))) < 1e-8)

## ---- The covariance structure is untouched by the mean structure: fitting
##      free means on raw data equals fitting stationary means on per-wave
##      mean-centered data ----
dcent <- s$data
for (t in seq_len(nT)) for (cc in s$design[, t]) dcent[[cc]] <- dcent[[cc]] - mean(dcent[[cc]])
mod_cent <- suppressWarnings(runmodel(panelgvar(dcent, vars = s$design), verbose = FALSE))
expect_equal(mod_free@fitmeasures$chisq, mod_cent@fitmeasures$chisq, tolerance = 1e-6)
for (mat in c("beta", "omega_zeta_within", "omega_zeta_between")){
  expect_true(max(abs(getmatrix(mod_free, mat) - getmatrix(mod_cent, mat))) < 1e-6)
}
# ... while stationary means on the raw (shifted) data misfit badly:
mod_raw <- suppressWarnings(runmodel(panelgvar(s$data, vars = s$design), verbose = FALSE))
expect_true(mod_raw@fitmeasures$chisq > 10 * mod_free@fitmeasures$chisq)

## ---- Baseline mirrors the target's mean structure; the saturated model
##      (which already has saturated means) is unchanged ----
expect_equal(max(mod_free@baseline_saturated$baseline@parameters$par),
             max(mod_stat@baseline_saturated$baseline@parameters$par) + m * (nT - 1))
expect_equal(max(mod_free@baseline_saturated$saturated@parameters$par),
             max(mod_stat@baseline_saturated$saturated@parameters$par))

## ---- Cross-framework equivalence: dlvm1 with lambda = I and zero residual
##      variances equals panelgvar, also with free means ----
O <- matrix(0, m, m)
mod_dl <- suppressWarnings(runmodel(dlvm1(s$data, vars = s$design, lambda = diag(m),
                                          within_latent = "ggm", between_latent = "ggm",
                                          sigma_epsilon_within = O, sigma_epsilon_between = O,
                                          omega_epsilon_within = O, omega_epsilon_between = O,
                                          delta_epsilon_within = O, delta_epsilon_between = O,
                                          stationary_means = FALSE), verbose = FALSE))
expect_equal(mod_dl@fitmeasures$df, mod_free@fitmeasures$df)
expect_equal(mod_dl@fitmeasures$npar, mod_free@fitmeasures$npar)
expect_true(abs(mod_dl@fitmeasures$logl - mod_free@fitmeasures$logl) < 1e-6)
expect_true(abs(mod_dl@fitmeasures$chisq - mod_free@fitmeasures$chisq) < 1e-6)
# mu_t(panelvar) = nu_t(dlvm1), with the latent means fixed to zero:
nu <- getmatrix(mod_dl, "nu"); mu_eta <- getmatrix(mod_dl, "mu_eta")
expect_equal(dim(nu), c(m, nT))
expect_true(all(mu_eta == 0))
expect_true(max(abs(matrix(getmatrix(mod_free, "mu"), m, nT) - nu)) < 1e-6)

## ---- dlvm1 with latent variables: wave-specific intercepts (m(T-1) extra
##      parameters), latent means still fixed to zero ----
set.seed(5)
q <- 2; nn <- 4; N <- 400
Lam <- matrix(0, nn, q); Lam[1:2, 1] <- c(1, .8); Lam[3:4, 2] <- c(1, .7)
Bq <- matrix(c(.3, .1, -.05, .25), q, q)
SigW <- matrix(c(.6, .15, .15, .5), q, q); SigZb <- matrix(c(.4, .1, .1, .3), q, q)
Sig0q <- matrix(solve(diag(q^2) - kronecker(Bq, Bq), as.vector(SigW)), q, q)
bi <- rmvnorm(N, rep(0, q), SigZb); ebi <- rmvnorm(N, rep(0, nn), diag(.2, nn))
eta <- rmvnorm(N, rep(0, q), Sig0q)
latshift <- c(0, .4, .8)
Yl <- matrix(0, N, nn * nT)
for (t in seq_len(nT)){
  if (t > 1) eta <- eta %*% t(Bq) + rmvnorm(N, rep(0, q), SigW)
  Yl[, (t - 1) * nn + seq_len(nn)] <- (eta + bi + latshift[t]) %*% t(Lam) + ebi +
    rmvnorm(N, rep(0, nn), diag(.4, nn))
}
colnames(Yl) <- paste0("y", rep(seq_len(nn), nT), "_", rep(seq_len(nT), each = nn))
dl <- as.data.frame(Yl)
designL <- matrix(colnames(dl), nn, nT); rownames(designL) <- paste0("y", seq_len(nn))
lampat <- 1 * (Lam != 0)

dl_stat <- suppressWarnings(runmodel(dlvm1(dl, vars = designL, lambda = lampat),
                                     verbose = FALSE))
dl_free <- suppressWarnings(runmodel(dlvm1(dl, vars = designL, lambda = lampat,
                                           stationary_means = FALSE), verbose = FALSE))
expect_equal(dl_free@fitmeasures$npar, dl_stat@fitmeasures$npar + nn * (nT - 1))
expect_equal(dl_stat@fitmeasures$df - dl_free@fitmeasures$df, nn * (nT - 1))
# The intercepts are wave-specific (one per indicator per wave) and all free ...
nup <- dl_free@parameters[dl_free@parameters$matrix == "nu", ]
expect_equal(nrow(nup), nn * nT)
expect_true(all(!nup$fixed))
# ... the latent means remain a single column, fixed to zero ...
me <- dl_free@parameters[dl_free@parameters$matrix == "mu_eta", ]
expect_equal(nrow(me), q)
expect_true(all(me$fixed))
expect_true(all(me$est == 0))
# ... and the mean structure is saturated, so the intercepts equal the observed
# per-wave means:
expect_true(max(abs(getmatrix(dl_free, "nu") - matrix(colMeans(dl), nn, nT))) < 1e-6)
# The implied mean differences recover the simulated level shifts (the shift
# enters indicator j scaled by its loading):
nu_hat <- getmatrix(dl_free, "nu")
lam_hat <- getmatrix(dl_free, "lambda")
shift_hat <- (nu_hat[, -1] - nu_hat[, 1]) / rowSums(lam_hat)
expect_true(max(abs(shift_hat - matrix(latshift[-1], nn, nT - 1, byrow = TRUE))) < 0.2)

## ---- The baseline model has the same (saturated, wave-specific) mean
##      structure as the latent model, so it is nested in it ----
base_free <- dl_free@baseline_saturated$baseline@parameters
expect_equal(sum(base_free$matrix == "mu" & !base_free$fixed), nn * nT)
expect_equal(sum(!nup$fixed), nn * nT)

## ---- The R and the C++ estimation paths agree (the C++ path is the default) ----
p_R <- suppressWarnings(runmodel(usecpp(panelgvar(s$data, vars = s$design, stationary_means = FALSE),
                                        FALSE), verbose = FALSE))
expect_true(abs(p_R@fitmeasures$logl - mod_free@fitmeasures$logl) < 1e-6)
expect_true(abs(p_R@fitmeasures$chisq - mod_free@fitmeasures$chisq) < 1e-6)
expect_true(max(abs(p_R@parameters$est - mod_free@parameters$est)) < 1e-6)

d_R <- suppressWarnings(runmodel(usecpp(dlvm1(dl, vars = designL, lambda = lampat,
                                              stationary_means = FALSE), FALSE), verbose = FALSE))
expect_true(abs(d_R@fitmeasures$logl - dl_free@fitmeasures$logl) < 1e-6)
expect_true(max(abs(d_R@parameters$est - dl_free@parameters$est)) < 1e-6)

## ---- Analytic and numeric gradients agree for the new mean part, in both
##      implementations (at_home: numDeriv on the full parameter vector) ----
if (at_home()){
  for (mod in list(panelgvar(s$data, vars = s$design, stationary_means = FALSE),
                   usecpp(panelgvar(s$data, vars = s$design, stationary_means = FALSE), FALSE),
                   dlvm1(dl, vars = designL, lambda = lampat, stationary_means = FALSE),
                   usecpp(dlvm1(dl, vars = designL, lambda = lampat, stationary_means = FALSE), FALSE))){
    cj <- suppressWarnings(psychonetrics:::checkJacobian(mod, plot = FALSE))
    expect_true(max(abs(cj$analytic - cj$numeric)) / max(abs(cj$numeric)) < 1e-4)
  }
}

## ---- Multiple groups: equal = "mu" equates the means per wave ----
s2 <- sim_panel(seed = 2, shifts = c(0, 0.5, 1))
dgrp <- rbind(cbind(s$data, group = "A"), cbind(s2$data, group = "B"))
g_free <- suppressWarnings(panelgvar(dgrp, vars = s$design, groupvar = "group",
                                     stationary_means = FALSE, equal = "mu"))
mu_g <- g_free@parameters[g_free@parameters$matrix == "mu", ]
expect_equal(nrow(mu_g), 2 * m * nT)
# Every (row, col) combination shares one parameter number across the groups:
expect_equal(length(unique(mu_g$par)), m * nT)
for (r in seq_len(m)) for (cc in seq_len(nT)){
  expect_equal(length(unique(mu_g$par[mu_g$row == r & mu_g$col == cc])), 1L)
}
