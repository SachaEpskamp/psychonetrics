# Tests for the composite variable model family (cvm / cnm), the factor- and
# composite-based SEM of Schamberger, Schuberth, Henseler & Rosseel.
#
# The tests avoid hard-coded reference values from simulated data: they check
# exact recovery from a population covariance matrix, agreement between the R
# and C++ pipelines, agreement between equivalent parameterizations, the
# analytic derivatives, and (when available) agreement with lavaan.

suppressMessages(library(psychonetrics))

# Shared helpers (portable multivariate normal draws; see _helpers.R):
source(if (file.exists("_helpers.R")) "_helpers.R" else
         system.file("tinytest", "_helpers.R", package = "psychonetrics"))

## ---- Population model (scenario of Schamberger et al.) ----
# Two common factors (f1: y1-y3, f2: y4-y6) and two composites (c3: y7-y9,
# c4: y10-y13); c4 ~ f1 + c3, f2 ~ f1 + c3 + c4, f1 ~~ c3.
Lf <- matrix(0, 6, 2); Lf[1:3, 1] <- c(1, .8, .9); Lf[4:6, 2] <- c(1, .7, .8)
Thf <- diag(.5, 6)
W <- matrix(0, 7, 2); W[1:3, 1] <- c(1, .4, .6); W[4:7, 2] <- c(1, .2, .5, .3)
T1 <- matrix(c(6, 2, 1.2, 2, 5, 1.5, 1.2, 1.5, 2), 3)
T2 <- matrix(c(7, .9, .5, 1.2, .9, 3, .7, 1.8, .5, .7, 2, .5, 1.2, 1.8, .5, 5), 4)
Tm <- matrix(0, 7, 7); Tm[1:3, 1:3] <- T1; Tm[4:7, 4:7] <- T2
B <- matrix(0, 4, 4); B[4, 1] <- .6; B[4, 3] <- .5; B[2, 1] <- .2; B[2, 3] <- .1; B[2, 4] <- .3
isC <- c(FALSE, FALSE, TRUE, TRUE)
D <- t(W) %*% Tm %*% W
A <- solve(diag(4) - B)
Psi <- matrix(0, 4, 4); Psi[1, 1] <- 2; Psi[2, 2] <- 1.5; Psi[1, 3] <- Psi[3, 1] <- .1
# The residual variances of the composites follow from their variances:
diag(Psi)[isC] <- solve(A[isC, isC]^2, diag(D) - diag(A %*% Psi %*% t(A))[isC])
V <- A %*% Psi %*% t(A)
Lc <- Tm %*% W %*% solve(D)
L <- rbind(cbind(Lf, matrix(0, 6, 2)), cbind(matrix(0, 7, 2), Lc))
Theta <- matrix(0, 13, 13); Theta[1:6, 1:6] <- Thf; Theta[7:13, 7:13] <- Tm - Lc %*% D %*% t(Lc)
Sigma <- L %*% V %*% t(L) + Theta
nm <- paste0("y", 1:13)
dimnames(Sigma) <- list(nm, nm)
mu <- seq(0.5, 6.5, by = 0.5)

# Model matrices:
lambda <- matrix(0, 13, 4); lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1
weights <- matrix(0, 13, 4); weights[7:9, 3] <- 1; weights[10:13, 4] <- 1
beta <- matrix(0, 4, 4); beta[4, 1] <- 1; beta[4, 3] <- 1; beta[2, 1] <- 1; beta[2, 3] <- 1; beta[2, 4] <- 1
sigma_zeta <- diag(4); sigma_zeta[1, 3] <- sigma_zeta[3, 1] <- 1
latents <- c("f1", "f2", "c3", "c4")

## ---- 1. Exact recovery from the population moments ----
mod_pop <- runmodel(cvm(covs = Sigma, means = mu, nobs = 1000, covtype = "ML",
                        lambda = lambda, weights = weights, beta = beta,
                        sigma_zeta = sigma_zeta, latents = latents))
expect_inherits(mod_pop, "psychonetrics")
expect_equal(mod_pop@fitmeasures$df, 52)
expect_true(abs(mod_pop@fitmeasures$chisq) < 1e-4)
expect_true(max(abs(getmatrix(mod_pop, "weights")[7:13, 3:4] - W)) < 1e-3)
expect_true(max(abs(getmatrix(mod_pop, "beta") - B)) < 1e-3)
expect_true(max(abs(getmatrix(mod_pop, "sigma") - Sigma)) < 1e-2)
# Computed (residual) variances of the composites:
expect_true(max(abs(diag(getmatrix(mod_pop, "sigma_zeta"))[isC] - diag(Psi)[isC])) < 1e-2)
# Composite loadings and the implied residual covariances of composite indicators:
expect_true(max(abs(getmatrix(mod_pop, "lambda")[7:13, 3:4] - Lc)) < 1e-3)
expect_true(max(abs(getmatrix(mod_pop, "sigma_epsilon_c") - Theta[7:13, 7:13])) < 1e-2)
# Implied means, and composite means equal to the weighted indicator means:
expect_true(max(abs(as.vector(getmatrix(mod_pop, "mu")) - mu)) < 1e-3)
expect_true(max(abs(as.vector(getmatrix(mod_pop, "mu_eta"))[isC] - as.vector(t(W) %*% mu[7:13]))) < 1e-2)

## ---- 2. Sample data: R and C++ pipelines agree ----
set.seed(1)
dat <- as.data.frame(rmvnorm(400, mu, Sigma))
names(dat) <- nm

mod <- cvm(dat, lambda = lambda, weights = weights, beta = beta,
           sigma_zeta = sigma_zeta, latents = latents)
fitC <- runmodel(mod)
fitR <- runmodel(usecpp(mod, FALSE))
expect_equal(fitC@fitmeasures$df, 52)
expect_equal(fitC@fitmeasures$npar, 52)
expect_true(max(abs(fitC@parameters$est - fitR@parameters$est)) < 1e-6)
expect_true(max(abs(fitC@parameters$se - fitR@parameters$se), na.rm = TRUE) < 1e-6)
expect_true(abs(fitC@fitmeasures$chisq - fitR@fitmeasures$chisq) < 1e-6)

# The variance of a composite equals w' T w:
Wh <- getmatrix(fitC, "weights")[7:13, 3:4]
expect_true(max(abs(diag(getmatrix(fitC, "sigma_eta"))[isC] -
                      diag(t(Wh) %*% getmatrix(fitC, "sigma_c") %*% Wh))) < 1e-8)

# Computed entries are fixed, and structural zeros can never be freed:
pars <- fitC@parameters
expect_true(all(pars$fixed[pars$matrix == "sigma_zeta" & pars$row == pars$col & pars$row %in% 3:4]))
expect_true(all(pars$fixed[pars$matrix == "nu_eta" & pars$row %in% 3:4]))
crossblock <- pars$matrix == "sigma_c" & (pars$var1 %in% nm[7:9]) != (pars$var2 %in% nm[7:9])
expect_true(all(pars$fixed[crossblock] & pars$identified[crossblock]))
expect_true(all(pars$identified[pars$matrix == "lambda" & pars$row > 6]))

## ---- 3. Analytic derivatives ----
x0 <- psychonetrics:::parVector(mod)
set.seed(2)
mod_pert <- psychonetrics:::updateModel(x0 + runif(length(x0), -0.05, 0.05), mod)
# (the R pipeline is only checked at home, to limit the run time on CRAN)
for (m in if (at_home()) list(mod_pert, usecpp(mod_pert, FALSE)) else list(mod_pert)){
  g <- checkJacobian(m, plot = FALSE)
  expect_true(max(abs(g$analytic - g$numeric)) < 1e-4)
  f <- checkFisher(m, plot = FALSE)
  expect_true(max(abs(f$analytic - f$numeric)) < 1e-4)
}
# The model Jacobian is the same in R and C++:
JR <- psychonetrics:::d_phi_theta_cvm(psychonetrics:::prepare_cvm(x0, usecpp(mod, FALSE)))
JC <- psychonetrics:::d_phi_theta_cvm_cpp(psychonetrics:::prepare_cvm_cpp(x0, mod))
expect_true(max(abs(as.matrix(JR) - as.matrix(JC))) < 1e-10)

## ---- 4. Equivalent parameterizations give the same fit ----
for (type in if (at_home()) c("chol", "prec", "ggm", "cor") else "chol"){
  fit_type <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                           sigma_zeta = sigma_zeta, latents = latents, composite = type))
  expect_equal(fit_type@fitmeasures$df, 52)
  expect_true(abs(fit_type@fitmeasures$chisq - fitC@fitmeasures$chisq) < 1e-3)
}

# Latent network (GGM) with the same structure: f1 -- c3 exogenous, f2 and c4 endogenous:
omega_zeta <- matrix(0, 4, 4); omega_zeta[1, 3] <- omega_zeta[3, 1] <- 1
fit_ggm <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                        omega_zeta = omega_zeta, latents = latents, latent = "ggm"))
expect_equal(fit_ggm@fitmeasures$df, 52)
expect_true(abs(fit_ggm@fitmeasures$chisq - fitC@fitmeasures$chisq) < 1e-3)
gg <- checkJacobian(psychonetrics:::updateModel(psychonetrics:::parVector(fit_ggm) * 1.02, fit_ggm), plot = FALSE)
expect_true(max(abs(gg$analytic - gg$numeric)) < 1e-4)

# A composite cannot have both incoming effects and edges in the latent network:
omega_full <- matrix(1, 4, 4); diag(omega_full) <- 0
expect_error(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                 omega_zeta = omega_full, latents = latents, latent = "ggm"))

# All constructs in one saturated latent network: omega_zeta holds the partial
# correlations of the construct covariances:
fit_net <- runmodel(cvm(dat, lambda = lambda, weights = weights, latents = latents, latent = "ggm"))
Vz <- getmatrix(fit_net, "sigma_zeta")
pc <- -cov2cor(solve(Vz)); diag(pc) <- 0
expect_true(max(abs(pc - getmatrix(fit_net, "omega_zeta"))) < 1e-8)

## ---- 5. Composite network model and model search ----
fit_cnm <- runmodel(cnm(dat, lambda = lambda, weights = weights, beta = beta,
                        sigma_zeta = sigma_zeta, latents = latents))
expect_equal(fit_cnm@submodel, "cnm")
expect_true(abs(fit_cnm@fitmeasures$chisq - fitC@fitmeasures$chisq) < 1e-3)
# omega_c equals the partial correlations of the composite indicator covariances:
Tc <- getmatrix(fit_cnm, "sigma_c")
pcT <- -cov2cor(solve(Tc)); diag(pcT) <- 0
expect_true(max(abs(pcT - getmatrix(fit_cnm, "omega_c"))) < 1e-8)
# No edges between indicators of different composites:
expect_true(all(getmatrix(fit_cnm, "omega_c")[1:3, 4:7] == 0))

fit_pruned <- suppressWarnings(stepup(prune(fit_cnm, alpha = 0.01)))
expect_true(fit_pruned@fitmeasures$df >= fit_cnm@fitmeasures$df)
expect_true(all(getmatrix(fit_pruned, "omega_c")[1:3, 4:7] == 0))
expect_error(cnm(dat, weights = weights, composite = "cov"))

## ---- 6. Input checks ----
bad <- weights; bad[1, 3] <- 1 # y1 is a factor indicator
expect_error(cvm(dat, lambda = lambda, weights = bad, latents = latents))
bad <- weights; bad[7, 4] <- 1 # y7 in two composites
expect_error(cvm(dat, lambda = lambda, weights = bad, latents = latents))
bad <- lambda; bad[1, 3] <- 1 # c3 with both loadings and weights
expect_error(cvm(dat, lambda = bad, weights = weights, latents = latents))

# An isolated composite with estimated weights gives a warning when the model is run:
w2 <- matrix(0, 7, 2); w2[1:3, 1] <- 1; w2[4:7, 2] <- 1
expect_warning(psychonetrics:::cvm_update_structure(cvm(dat[, 7:13], weights = w2, sigma_zeta = "diag"), warn = TRUE))

## ---- 7. Without composites the model equals an lvm ----
lam2 <- matrix(0, 6, 2); lam2[1:3, 1] <- 1; lam2[4:6, 2] <- 1
fit_c <- runmodel(cvm(dat[, 1:6], lambda = lam2))
fit_l <- runmodel(lvm(dat[, 1:6], lambda = lam2))
expect_equal(fit_c@fitmeasures$df, fit_l@fitmeasures$df)
expect_true(abs(fit_c@fitmeasures$chisq - fit_l@fitmeasures$chisq) < 1e-2)

## ---- 8. Scaling: a first indicator with a (near) zero weight ----
# Composite cc of y4-y6 with weights (0, 1, .6): the first indicator is a bad
# choice to scale the composite, so another indicator must be used instead.
Tb <- matrix(c(2, .6, .4, .6, 1.5, .5, .4, .5, 1), 3)
wb <- c(0, 1, .6)
db <- drop(t(wb) %*% Tb %*% wb)
bb <- 0.5 / sqrt(db)
Ab <- matrix(c(1, 0, bb, 1), 2) # constructs (f1, cc): f1 <- cc
Vb <- Ab %*% diag(c(1 - bb^2 * db, db)) %*% t(Ab)
Lb <- rbind(cbind(c(1, .8, .9), 0), cbind(0, Tb %*% wb / db))
Thb <- matrix(0, 6, 6); diag(Thb)[1:3] <- .5
Thb[4:6, 4:6] <- Tb - (Tb %*% wb) %*% t(Tb %*% wb) / db
Sb <- Lb %*% Vb %*% t(Lb) + Thb
dimnames(Sb) <- list(nm[1:6], nm[1:6])
lamb <- matrix(0, 6, 2); lamb[1:3, 1] <- 1
wtsb <- matrix(0, 6, 2); wtsb[4:6, 2] <- 1
betab <- matrix(0, 2, 2); betab[1, 2] <- 1

expect_message(mod_bm <- cvm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                             beta = betab, latents = c("f1", "cc")), "scale")
pb <- mod_bm@parameters
# The second indicator (y5) scales the composite and the weight of the first is free:
expect_equal(pb$row[pb$matrix == "weights" & pb$col == 2 & pb$fixed & pb$est == 1], 5)
expect_true(!pb$fixed[pb$matrix == "weights" & pb$col == 2 & pb$row == 4])
fit_bm <- runmodel(mod_bm)
expect_true(abs(fit_bm@fitmeasures$chisq) < 1e-4)
expect_true(max(abs(getmatrix(fit_bm, "weights")[4:6, 2] - wb)) < 1e-2)
expect_equal(fit_bm@optim$convergence, 0)

# The first indicator is always used with bad_marker_crit = 0:
mod_first <- cvm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                 beta = betab, latents = c("f1", "cc"), bad_marker_crit = 0)
pf <- mod_first@parameters
expect_equal(pf$row[pf$matrix == "weights" & pf$col == 2 & pf$fixed & pf$est == 1], 4)

# The optimizer converges for the models above:
expect_equal(fitC@optim$convergence, 0)
expect_equal(fit_ggm@optim$convergence, 0)

## ---- 9. Longer-running checks ----
if (at_home()){

  # An observed variable in the structural model (single-indicator composite),
  # multiple groups with equality constraints, and missing data:
  dat$group <- rep(c("a", "b"), each = 200)
  fit_conf <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta,
                           latents = latents, vars = nm, groupvar = "group"))
  fit_weak <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta,
                           latents = latents, vars = nm, groupvar = "group", equal = c("lambda", "weights")))
  fit_strong <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta,
                             latents = latents, vars = nm, groupvar = "group", equal = c("lambda", "weights", "nu")))
  expect_equal(fit_conf@fitmeasures$df, 104)
  expect_equal(fit_weak@fitmeasures$df, 113)
  expect_equal(fit_strong@fitmeasures$df, 124)
  expect_true(fit_weak@fitmeasures$chisq >= fit_conf@fitmeasures$chisq - 1e-4)
  expect_true(fit_strong@fitmeasures$chisq >= fit_weak@fitmeasures$chisq - 1e-4)
  fit_weakR <- runmodel(usecpp(cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta,
                                   latents = latents, vars = nm, groupvar = "group", equal = c("lambda", "weights")), FALSE))
  expect_true(max(abs(fit_weak@parameters$est - fit_weakR@parameters$est)) < 1e-6)

  dat_mis <- dat[, nm]
  set.seed(3)
  for (j in c(2, 8, 11)) dat_mis[sample(400, 40), j] <- NA
  mod_fiml <- cvm(dat_mis, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta, latents = latents)
  expect_equal(mod_fiml@estimator, "FIML")
  fit_fiml <- runmodel(mod_fiml)
  fit_fimlR <- runmodel(usecpp(mod_fiml, FALSE))
  expect_equal(fit_fiml@fitmeasures$df, 52)
  expect_true(max(abs(fit_fiml@parameters$est - fit_fimlR@parameters$est)) < 1e-6)

  # Comparison with lavaan (composites via the "<~" operator):
  if (requireNamespace("lavaan", quietly = TRUE)){
    lav_model <- ' f1 =~ y1 + y2 + y3
                   f2 =~ y4 + y5 + y6
                   c3 <~ y7 + y8 + y9
                   c4 <~ y10 + y11 + y12 + y13
                   c4 ~ f1 + c3
                   f2 ~ f1 + c3 + c4
                   f1 ~~ c3 '
    lav_fit <- tryCatch(lavaan::sem(lav_model, data = dat[, nm], meanstructure = TRUE, likelihood = "normal",
                                    composites.cov = "free"),
                        error = function(e) NULL)
    if (!is.null(lav_fit) && isTRUE(lavaan::lavInspect(lav_fit, "converged"))){
      expect_equal(as.numeric(lavaan::fitMeasures(lav_fit, "df")), fitC@fitmeasures$df)
      expect_true(abs(as.numeric(lavaan::fitMeasures(lav_fit, "chisq")) - fitC@fitmeasures$chisq) < 1e-3)
      expect_true(abs(as.numeric(lavaan::fitMeasures(lav_fit, "logl")) - fitC@fitmeasures$logl) < 1e-3)
      pe <- lavaan::parameterEstimates(lav_fit)
      w_psy <- fitC@parameters[fitC@parameters$matrix == "weights" & fitC@parameters$est != 0, ]
      expect_true(max(abs(w_psy$est - pe$est[pe$op == "<~"])) < 1e-2)
      expect_true(max(abs(w_psy$se - pe$se[pe$op == "<~"]), na.rm = TRUE) < 1e-2)
      expect_true(max(abs(getmatrix(fitC, "sigma") - lavaan::lavInspect(lav_fit, "implied")$cov[nm, nm])) < 1e-2)
    }
  }
}
