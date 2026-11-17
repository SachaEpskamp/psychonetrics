# Tests for the composite variable model family (cvm / cnm), the factor- and
# composite-based SEM of Schamberger, Schuberth, Henseler & Rosseel.
#
# The tests avoid hard-coded reference values from simulated data: they check
# exact recovery from a population covariance matrix, agreement between the R
# and C++ pipelines, agreement between equivalent parameterizations, the
# analytic derivatives, and (when available) agreement with lavaan.
#
# Tolerances: the fit function of a composite model is flat in the direction of
# the composite weights, and the package uses the same optimizer settings for
# all model families. Fits therefore stop slightly short of the optimum
# (typically 1e-4 to 1e-3 in the chi-square here). Comparisons of separately
# optimized fits use tolerances with a wide margin around that; comparisons
# that do not depend on the optimizer (R versus C++ at the same point,
# derivatives) remain strict.

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
expect_true(abs(mod_pop@fitmeasures$chisq) < 1e-2)
expect_true(max(abs(getmatrix(mod_pop, "weights")[7:13, 3:4] - W)) < 1e-2)
expect_true(max(abs(getmatrix(mod_pop, "beta") - B)) < 1e-2)
expect_true(max(abs(getmatrix(mod_pop, "sigma") - Sigma)) < 5e-2)
# Computed (residual) variances of the composites:
expect_true(max(abs(diag(getmatrix(mod_pop, "sigma_zeta"))[isC] - diag(Psi)[isC])) < 1e-1)
# Composite loadings and the implied residual covariances of composite indicators:
expect_true(max(abs(getmatrix(mod_pop, "lambda")[7:13, 3:4] - Lc)) < 1e-2)
expect_true(max(abs(getmatrix(mod_pop, "sigma_epsilon_c") - Theta[7:13, 7:13])) < 1e-2)
# Implied means, and composite means equal to the weighted indicator means:
expect_true(max(abs(as.vector(getmatrix(mod_pop, "mu")) - mu)) < 1e-2)
expect_true(max(abs(as.vector(getmatrix(mod_pop, "mu_eta"))[isC] - as.vector(t(W) %*% mu[7:13]))) < 5e-2)

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
  expect_true(abs(fit_type@fitmeasures$chisq - fitC@fitmeasures$chisq) < 5e-2)
}

# Latent network (GGM) with the same structure: f1 -- c3 exogenous, f2 and c4 endogenous:
omega_zeta <- matrix(0, 4, 4); omega_zeta[1, 3] <- omega_zeta[3, 1] <- 1
fit_ggm <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                        omega_zeta = omega_zeta, latents = latents, latent = "ggm"))
expect_equal(fit_ggm@fitmeasures$df, 52)
expect_true(abs(fit_ggm@fitmeasures$chisq - fitC@fitmeasures$chisq) < 5e-2)
gg <- checkJacobian(psychonetrics:::updateModel(psychonetrics:::parVector(fit_ggm) * 1.02, fit_ggm), plot = FALSE)
expect_true(max(abs(gg$analytic - gg$numeric)) < 1e-4)

# A composite cannot have both incoming effects and edges in the latent network:
omega_full <- matrix(1, 4, 4); diag(omega_full) <- 0
expect_error(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                 omega_zeta = omega_full, latents = latents, latent = "ggm"))

# The identification step marks the parameters that would give such a model, so
# that model search skips them. In fit_ggm, c3 is a node of the latent network
# (no incoming effects) and c4 has incoming effects: effects on c3 and edges of
# c4 are marked; the common factors are not restricted:
Pg <- fit_ggm@parameters
expect_true(all(Pg$identified[Pg$matrix == "beta" & Pg$fixed & Pg$var1 == "c3" & Pg$var2 != "c3"]))
expect_true(all(Pg$identified[Pg$matrix == "omega_zeta" & Pg$fixed & Pg$var1 != Pg$var2 & (Pg$var1 == "c4" | Pg$var2 == "c4")]))
expect_false(any(Pg$identified[Pg$matrix == "omega_zeta" & Pg$fixed & Pg$var1 == "f2" & Pg$var2 == "f1"]))
expect_false(any(Pg$identified[Pg$matrix == "beta" & Pg$fixed & Pg$var1 == "f1" & Pg$var2 != "f1"]))

if (at_home()){
  # Model search over effects and edges together does not give such a model
  # (runmodel() itself contains no check: the marks come from identification):
  beta_s <- matrix(0, 4, 4); beta_s[4, 3] <- 1; beta_s[2, 1] <- 1
  fit_s <- suppressWarnings(stepup(runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta_s,
                                                omega_zeta = "zero", latents = latents, latent = "ggm")),
                                   matrices = c("beta", "omega_zeta"), alpha = 0.01))
  Bs <- getmatrix(fit_s, "beta"); Os <- getmatrix(fit_s, "omega_zeta")
  expect_false(any(rowSums(Bs[3:4, , drop = FALSE] != 0) > 0 & rowSums(Os[3:4, , drop = FALSE] != 0) > 0))
  expect_true(fit_s@fitmeasures$df < 57)

  # An effect that is freed by hand on a composite that is a node of the network:
  # the implied model refuses it when the model is run, and identify() explains why:
  bad_ggm <- freepar(fit_ggm, "beta", "c3", "f2")
  expect_error(suppressWarnings(runmodel(bad_ggm)), "both incoming effects")
  expect_error(identify(bad_ggm), "both incoming effects")

  # The computed variance of a composite is fixed by the identification step. Freed
  # by hand it gives a warning after the model is run; identify() fixes it again:
  free_var <- freepar(fitC, "sigma_zeta", "c3", "c3")
  expect_warning(runmodel(free_var), "cannot be estimated")
  Pf <- identify(free_var)@parameters
  expect_true(all(Pf$fixed[Pf$matrix == "sigma_zeta" & Pf$var1 == "c3" & Pf$var2 == "c3"]))
}

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
expect_true(abs(fit_cnm@fitmeasures$chisq - fitC@fitmeasures$chisq) < 5e-2)
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

## ---- 5b. transmod(): switching between equivalent cvm and cnm models ----
# The '_c' matrices are block diagonal (one block per composite). Every
# parameterization of a block diagonal covariance matrix is block diagonal, so
# a model with saturated blocks can be transformed between the types; the
# elements between indicators of different composites stay fixed to zero.
tab_structure <- function(m){
  P <- m@parameters
  sort(paste(P$matrix, P$var1, P$var2, P$var1_id, P$var2_id, P$group_id, P$fixed, P$identified, P$par > 0))
}
objective_at <- function(m) psychonetrics:::psychonetrics_fitfunction_cpp(psychonetrics:::parVector(m), m)

# cvm -> cnm: the same model as specified with cnm(), at the same solution:
to_cnm <- transmod(fitC, composite = "ggm")
expect_equal(to_cnm@types$composite, "ggm")
expect_equal(to_cnm@submodel, "cnm")
expect_identical(tab_structure(to_cnm), tab_structure(fit_cnm))
expect_true(abs(objective_at(to_cnm) - fitC@objective) < 1e-8)
pc0 <- -cov2cor(solve(getmatrix(fitC, "sigma_c"))); diag(pc0) <- 0
expect_true(max(abs(pc0 - suppressWarnings(getmatrix(to_cnm, "omega_c")))) < 1e-8)

# cnm -> cvm:
to_cvm <- transmod(fit_cnm, composite = "cov")
expect_equal(to_cvm@submodel, "cvm")
expect_identical(tab_structure(to_cvm), tab_structure(fitC))
expect_true(abs(objective_at(to_cvm) - fit_cnm@objective) < 1e-8)

# The residual type (a diagonal structure): no parameters are added, same fit:
to_res <- transmod(fitC, residual = "ggm")
expect_equal(max(to_res@parameters$par), max(fitC@parameters$par))
expect_false(any(to_res@parameters$fixed & to_res@parameters$par > 0))
expect_true(abs(objective_at(to_res) - fitC@objective) < 1e-8)

# Not possible: other latent types than cov / ggm, and structures that are
# neither saturated nor diagonal:
expect_error(transmod(fitC, latent = "chol"), "can only be 'cov' or 'ggm'")
expect_error(transmod(fit_pruned, composite = "cov"), "Only diagonal and saturated")

if (at_home()){
  for (type in c("chol", "prec", "cor")){
    expect_true(abs(objective_at(transmod(fitC, composite = type)) - fitC@objective) < 1e-8)
  }
  # Running the transformed model gives the same fit:
  run_cnm <- runmodel(to_cnm)
  expect_equal(run_cnm@fitmeasures$df, fitC@fitmeasures$df)
  expect_true(abs(run_cnm@fitmeasures$chisq - fitC@fitmeasures$chisq) < 5e-2)
  # With keep_computed = TRUE the standard errors are those of the cnm model:
  kept <- transmod(fitC, composite = "ggm", keep_computed = TRUE)
  sel_k <- kept@parameters$matrix == "omega_c" & !kept@parameters$fixed
  expect_true(all(is.finite(kept@parameters$se[sel_k])))
  expect_true(max(abs(kept@parameters$se[sel_k] - fit_cnm@parameters$se[sel_k])) < 1e-2)
  # A saturated latent network to covariances between the constructs: same
  # fit, and the marks of the network restriction are removed:
  to_cov <- transmod(fit_net, latent = "cov")
  expect_true(abs(objective_at(to_cov) - fit_net@objective) < 1e-8)
  expect_false(any(to_cov@parameters$identified[to_cov@parameters$matrix == "beta" &
                                                  to_cov@parameters$row != to_cov@parameters$col]))
}

## ---- 6. Input checks ----
bad <- weights; bad[1, 3] <- 1 # y1 is a factor indicator
expect_error(cvm(dat, lambda = lambda, weights = bad, latents = latents))
bad <- weights; bad[7, 4] <- 1 # y7 in two composites
expect_error(cvm(dat, lambda = lambda, weights = bad, latents = latents))
bad <- lambda; bad[1, 3] <- 1 # c3 with both loadings and weights
expect_error(cvm(dat, lambda = bad, weights = weights, latents = latents), "single indicator")

# An isolated composite with estimated weights gives a warning after the model is run
# (one of the warnings on the solution):
w2 <- matrix(0, 7, 2); w2[1:3, 1] <- 1; w2[4:7, 2] <- 1
expect_warning(psychonetrics:::cvm_solution_warnings(cvm(dat[, 7:13], weights = w2, sigma_zeta = "diag")),
               "not related to any other construct")

## ---- 6b. Starting values of the intercepts of factor indicators ----
# f2 is predicted by composites, which have nonzero means, so that f2 has a
# nonzero implied mean and the intercepts of its indicators (y4-y6) are not
# their means. The starting values take this into account: the implied means
# at the starting values equal the sample means. (With intercepts started at
# the sample means, ULS and DWLS fits of such models did not converge.)
mod_start <- cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta, latents = latents)
mu_start <- as.vector(as.matrix(psychonetrics:::implied_cvm(mod_start, all = FALSE)[[1]]$mu))
expect_true(max(abs(mu_start - colMeans(dat[, nm]))) < 1e-8)
nu_start <- mod_start@parameters$est[mod_start@parameters$matrix == "nu"]
expect_true(max(abs(nu_start[4:6] - colMeans(dat[, nm[4:6]]))) > 0.1)
expect_true(max(abs(nu_start[c(1:3, 7:13)] - colMeans(dat[, nm[c(1:3, 7:13)]]))) < 1e-8)
for (est in if (at_home()) c("DWLS", "ULS") else "DWLS"){
  fit_ls <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta, sigma_zeta = sigma_zeta,
                         latents = latents, estimator = est))
  expect_equal(fit_ls@optim$convergence, 0)
  # A converged solution: no gradient left, and reproducible from the R pipeline:
  expect_true(max(abs(psychonetrics:::psychonetrics_gradient_cpp(psychonetrics:::parVector(fit_ls), fit_ls))) < 1e-2)
  expect_true(all(diag(getmatrix(fit_ls, "sigma_zeta")) > 0))
}

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

# By default (bad_marker_crit = 1) a composite is scaled by the indicator with
# the largest standardized starting weight, and one message names the
# composites for which that is not the first indicator. With a smaller value
# the first indicator is only replaced when its starting weight is below that
# proportion of the largest, and the message states that rule:
expect_message(cvm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                   beta = betab, latents = c("f1", "cc")), "largest \\(standardized\\) starting weight to 1.*cc \\(y5\\)")
expect_message(cvm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                   beta = betab, latents = c("f1", "cc"), bad_marker_crit = 0.1), "smaller than 0.1 times")

# The family is flagged experimental with the shared note, which is shown once
# per session and is silent from version 0.19 on:
if (utils::packageVersion("psychonetrics") < "0.19"){
  exp_env <- psychonetrics:::.experimental_warned
  rm(list = grep("cvm", ls(exp_env), value = TRUE), envir = exp_env)
  expect_message(cnm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                     beta = betab, latents = c("f1", "cc")), "is experimental in psychonetrics")
  msgs_exp <- character(0)
  withCallingHandlers(
    cvm(covs = Sb, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
        beta = betab, latents = c("f1", "cc")),
    message = function(m){ msgs_exp <<- c(msgs_exp, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_false(any(grepl("experimental", msgs_exp)))
}
# Composite c3 of the first model (weights 1, .4, .6 with standard deviations
# 2.45, 2.24 and 1.41) has its largest standardized weight on the first
# indicator, which is then kept:
expect_equal(pars$row[pars$matrix == "weights" & pars$col == 3 & pars$fixed & pars$est == 1], 7)

## ---- 8a. Scaling: a composite that is reversed relative to its scaling indicator ----
# Composite cc of y4-y6 with weights (-.1, 1, .6): the first indicator has a
# small negative weight but correlates positively with the composite. Scaled
# by that indicator (weight fixed to 1) the composite is reversed, which gives
# a warning after the model is run; the default scaling does not:
wr <- c(-.1, 1, .6)
dr <- drop(t(wr) %*% Tb %*% wr)
br <- 0.5 / sqrt(dr)
Ar <- matrix(c(1, 0, br, 1), 2)
Vr <- Ar %*% diag(c(1 - br^2 * dr, dr)) %*% t(Ar)
Lr <- rbind(cbind(c(1, .8, .9), 0), cbind(0, Tb %*% wr / dr))
Thr <- matrix(0, 6, 6); diag(Thr)[1:3] <- .5
Thr[4:6, 4:6] <- Tb - (Tb %*% wr) %*% t(Tb %*% wr) / dr
Sr <- Lr %*% Vr %*% t(Lr) + Thr
dimnames(Sr) <- list(nm[1:6], nm[1:6])
expect_warning(fit_rev <- runmodel(cvm(covs = Sr, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                                       beta = betab, latents = c("f1", "cc"), bad_marker_crit = 0)),
               "correlates negatively")
expect_true(getmatrix(fit_rev, "lambda")[4, 2] < 0)
fit_notrev <- suppressMessages(runmodel(cvm(covs = Sr, nobs = 1000, covtype = "ML", lambda = lamb, weights = wtsb,
                                            beta = betab, latents = c("f1", "cc"))))
expect_silent(psychonetrics:::cvm_solution_warnings(fit_notrev))
expect_true(abs(fit_rev@fitmeasures$chisq - fit_notrev@fitmeasures$chisq) < 1e-2)

# A warning is also given when the scaling indicator has a small weight
# relative to the other indicators (a standardized weight below 1/25 of the
# largest). Here the composite of section 8 is scaled by y4 with weights
# (1, 50, 30), without estimating the model again:
mod_small <- fixpar(freepar(fit_bm, "weights", "y5", "cc", start = 50), "weights", "y4", "cc", value = 1)
mod_small@parameters$est[mod_small@parameters$matrix == "weights" & mod_small@parameters$var1 == "y6" &
                           mod_small@parameters$var2 == "cc"] <- 30
mod_small@modelmatrices <- psychonetrics:::formModelMatrices(mod_small)
expect_warning(psychonetrics:::cvm_solution_warnings(mod_small), "small weight relative to the other indicators")

# The optimizer converges for the models above:
expect_equal(fitC@optim$convergence, 0)
expect_equal(fit_ggm@optim$convergence, 0)

## ---- 8b. Penalized maximum likelihood (PML) ----
# By default the network of the composite indicators is penalized; edges are
# shrunk to exactly zero, the weights are never penalized, and refit() gives
# the unpenalized estimates of the selected network. (A small model with the
# two composites only: the proximal gradient optimizer used for PML can need
# very many iterations in larger composite models.)
weights_c <- matrix(0, 7, 2); weights_c[1:3, 1] <- 1; weights_c[4:7, 2] <- 1
beta_c <- matrix(0, 2, 2); beta_c[2, 1] <- 1
pml_model <- function(data = dat, ...){
  cnm(data, vars = nm[7:13], weights = weights_c, beta = beta_c, latents = c("c3", "c4"), ...)
}
mod_pml <- pml_model(estimator = "PML", penalty_lambda = 0.3)
expect_equal(mod_pml@estimator, "PML")
Ppml <- mod_pml@parameters
expect_true(all(Ppml$penalty_lambda[Ppml$matrix == "omega_c" & !Ppml$fixed] == 0.3))
expect_true(all(Ppml$penalty_lambda[Ppml$matrix != "omega_c"] == 0))
# The baseline and saturated models are not penalized:
expect_equal(mod_pml@baseline_saturated$saturated@estimator, "ML")
fit_pml <- suppressWarnings(suppressMessages(runmodel(mod_pml)))
expect_equal(fit_pml@optim$convergence, 0)
edges_pml <- fit_pml@parameters[fit_pml@parameters$matrix == "omega_c" & !fit_pml@parameters$fixed, ]
expect_true(any(edges_pml$est == 0))
expect_true(any(edges_pml$est != 0))
fit_refit <- suppressWarnings(suppressMessages(refit(fit_pml)))
fit_ml_c <- runmodel(pml_model())
expect_equal(fit_refit@estimator, "ML")
expect_equal(fit_refit@fitmeasures$df, fit_ml_c@fitmeasures$df + sum(edges_pml$est == 0))
expect_true(fit_refit@fitmeasures$chisq > fit_ml_c@fitmeasures$chisq - 1e-6)
# Without a network in the model there is nothing to penalize by default:
expect_warning(cvm(dat, vars = nm[7:13], weights = weights_c, beta = beta_c, latents = c("c3", "c4"),
                   estimator = "PML", penalty_lambda = 0.3), "No matrices to penalize")

## ---- 8c. A saturated construct part: closed-form starting values ----
# Three composites of three indicators in a saturated construct model (all
# composites related, no effects). With saturated covariances within the
# composites the chi-square is N (sum_k log|S_kk| + log|R_C| - log|S|), with
# R_C the correlation matrix of the composites: the maximum likelihood weights
# minimize |R_C| (GENVAR). The starting values are at that optimum:
Tp <- kronecker(diag(3), matrix(c(1, .4, .3, .4, 1, .2, .3, .2, 1), 3))
Wp <- matrix(0, 9, 3); Wp[1:3, 1] <- c(1, .8, .6); Wp[4:6, 2] <- c(1, .5, .7); Wp[7:9, 3] <- c(1, .6, .4)
Dp <- diag(diag(t(Wp) %*% Tp %*% Wp))
Rp <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3)
Lp <- Tp %*% Wp %*% solve(Dp)
Sp <- Tp + Lp %*% (sqrt(Dp) %*% Rp %*% sqrt(Dp) - Dp) %*% t(Lp)
set.seed(4)
dat3 <- as.data.frame(rmvnorm(400, rep(0, 9), Sp))
names(dat3) <- nm[1:9]
S3 <- cov(dat3) * 399 / 400
block3 <- rep(1:3, each = 3)
logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
logdet_RC <- function(w){
  Wm <- matrix(0, 9, 3); Wm[cbind(1:9, block3)] <- w
  logdet(cov2cor(t(Wm) %*% S3 %*% Wm))
}
opt3 <- optim(rep(1, 9), logdet_RC, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
chisq_closed <- 400 * (sum(sapply(1:3, function(k) logdet(S3[block3 == k, block3 == k]))) + opt3$value - logdet(S3))
mod3 <- suppressMessages(cvm(dat3, weights = 1 * (Wp != 0), latents = c("A", "B", "C")))
chisq_start <- 400 * (objective_at(mod3) - mod3@baseline_saturated$saturated@objective)
expect_true(abs(chisq_start - chisq_closed) < 1e-3)
fit3 <- runmodel(mod3)
expect_equal(fit3@fitmeasures$df, 18)
expect_true(abs(fit3@fitmeasures$chisq - chisq_closed) < 1e-3)
# The mode B weights that are used otherwise do not minimize the criterion:
adj3 <- matrix(TRUE, 3, 3); diag(adj3) <- FALSE
V_modeB <- psychonetrics:::cvm_start_weights(S3, block3, rep(TRUE, 3), matrix(0, 9, 3), adj3)
V_genvar <- psychonetrics:::cvm_start_weights(S3, block3, rep(TRUE, 3), matrix(0, 9, 3), adj3, genvar = TRUE)
expect_true(logdet(t(V_genvar) %*% S3 %*% V_genvar) < logdet(t(V_modeB) %*% S3 %*% V_modeB))
expect_true(abs(logdet(t(V_genvar) %*% S3 %*% V_genvar) - opt3$value) < 1e-6)

## ---- 8d. prune() can isolate a composite, with a warning ----
# The third composite is only weakly related to the other two (population
# covariances as input, so that its two edges are not significant). prune()
# removes both edges; the weights of the composite are then not identified,
# which is reported once:
Rw <- matrix(c(1, .5, .03, .5, 1, .02, .03, .02, 1), 3)
Sw <- Tp + Lp %*% (sqrt(Dp) %*% Rw %*% sqrt(Dp) - Dp) %*% t(Lp)
dimnames(Sw) <- list(nm[1:9], nm[1:9])
fit_w <- suppressMessages(suppressWarnings(runmodel(cnm(covs = Sw, nobs = 500, covtype = "ML", weights = 1 * (Wp != 0),
                                                        latent = "ggm", latents = c("A", "B", "C")))))
expect_warning(fit_iso <- prune(fit_w, alpha = 0.01, matrices = "omega_zeta"), "removed all relations.*: C\\.")
expect_true(all(getmatrix(fit_iso, "omega_zeta")[3, ] == 0))
expect_true(getmatrix(fit_iso, "omega_zeta")[2, 1] != 0)
expect_equal(psychonetrics:::cvm_unconnected_composites(fit_iso), "C")
if (at_home()){
  # With recursive pruning the warning is given once, and not again when a
  # model with an isolated composite is pruned further:
  warn_iso <- character(0)
  fit_iso2 <- withCallingHandlers(prune(fit_w, alpha = 0.01, matrices = "omega_zeta", recursive = TRUE),
                                  warning = function(w){ warn_iso <<- c(warn_iso, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_equal(sum(grepl("removed all relations", warn_iso)), 1)
  warn_iso <- character(0)
  fit_iso3 <- withCallingHandlers(prune(fit_iso, alpha = 0.01, matrices = "omega_c"),
                                  warning = function(w){ warn_iso <<- c(warn_iso, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_equal(sum(grepl("removed all relations", warn_iso)), 0)
}

## ---- 8e. Sum scores: all weights fixed to 1 (sumscores = TRUE) ----
# Composites c3 (y7-y9) and c4 (y10-y13) as sum scores that covary. Maximum
# likelihood then has a closed form: the covariances within the composites are
# the sample covariances, the covariances of the composites are those of the
# sum scores, and the chi-square is N * (sum_k log|S_kk| + log|R_C| - log|S|),
# with R_C the correlation matrix of the sum scores.
ss_vars <- nm[7:13]
ss_block <- rep(1:2, c(3, 4))
ss_w <- 1 * outer(ss_block, 1:2, "==")
ss_S <- cov(dat[, ss_vars]) * (nrow(dat) - 1) / nrow(dat)
msgs_ss <- character(0)
mod_ss <- withCallingHandlers(
  cvm(dat, vars = ss_vars, weights = ss_w, latents = c("c3", "c4"), sumscores = TRUE),
  message = function(m){ msgs_ss <<- c(msgs_ss, conditionMessage(m)); invokeRestart("muffleMessage") })
# No message on the scaling of the composites:
expect_false(any(grepl("scaled", msgs_ss)))
# All weights are fixed to 1, as restrictions (not as identification constraints):
p_ss <- mod_ss@parameters
w_ss <- p_ss$matrix == "weights" & p_ss$est != 0
expect_equal(sum(w_ss), 7)
expect_true(all(p_ss$fixed[w_ss] & p_ss$est[w_ss] == 1 & !p_ss$identified[w_ss]))
expect_false(any(p_ss$matrix == "weights" & !p_ss$fixed))
fit_ss <- runmodel(mod_ss)
expect_equal(fit_ss@fitmeasures$df, 7 * 8 / 2 - (3 * 4 / 2 + 4 * 5 / 2) - 1)
ss_chisq <- nrow(dat) * (sum(sapply(1:2, function(k) determinant(ss_S[ss_block == k, ss_block == k])$modulus)) +
                           determinant(cov2cor(t(ss_w) %*% ss_S %*% ss_w))$modulus - determinant(ss_S)$modulus)
expect_true(abs(fit_ss@fitmeasures$chisq - ss_chisq) < 1e-4)
expect_true(max(abs(getmatrix(fit_ss, "sigma_zeta") - t(ss_w) %*% ss_S %*% ss_w)) < 1e-4)
expect_true(all(as.matrix(getmatrix(fit_ss, "weights")) == ss_w))
expect_error(cvm(dat, vars = ss_vars, weights = ss_w, sumscores = "yes"), "sumscores")

# With latent = "ggm" the latent network is the network of the sum scores: the
# same partial correlations as ggm() on the sum scores, and the same chi-square
# difference when an edge is removed:
all_w <- 1 * outer(rep(1:4, c(3, 3, 3, 4)), 1:4, "==")
sumdat <- as.data.frame(as.matrix(dat[, nm]) %*% all_w)
names(sumdat) <- latents
om_ss <- 1 - diag(4); om_ss[4, 1] <- om_ss[1, 4] <- 0
fit_ssn <- runmodel(cvm(dat, vars = nm, weights = all_w, latents = latents, sumscores = TRUE, latent = "ggm"))
fit_ssn0 <- runmodel(cvm(dat, vars = nm, weights = all_w, latents = latents, sumscores = TRUE, latent = "ggm",
                         omega_zeta = om_ss))
fit_sumggm <- runmodel(ggm(sumdat))
fit_sumggm0 <- runmodel(ggm(sumdat, omega = om_ss))
expect_true(max(abs(getmatrix(fit_ssn, "omega_zeta") - getmatrix(fit_sumggm, "omega"))) < 1e-4)
expect_true(max(abs(getmatrix(fit_ssn0, "omega_zeta") - getmatrix(fit_sumggm0, "omega"))) < 1e-3)
expect_equal(fit_ssn0@fitmeasures$df, fit_ssn@fitmeasures$df + 1)
expect_true(abs((fit_ssn0@fitmeasures$chisq - fit_ssn@fitmeasures$chisq) - fit_sumggm0@fitmeasures$chisq) < 1e-3)

# With common factors and effects (the model of section 2): nested in the model
# with estimated weights (five weights fewer). Modification indices are
# computed for the fixed weights, and a weight can be freed:
fit_ss2 <- runmodel(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                        sigma_zeta = sigma_zeta, latents = latents, sumscores = TRUE))
expect_equal(fit_ss2@optim$convergence, 0)
expect_equal(fit_ss2@fitmeasures$df, fitC@fitmeasures$df + 5)
expect_true(fit_ss2@fitmeasures$chisq > fitC@fitmeasures$chisq)
p_ss2 <- fit_ss2@parameters
expect_true(all(is.finite(p_ss2$mi[p_ss2$matrix == "weights" & p_ss2$est != 0])))
fit_ss3 <- runmodel(freepar(fit_ss2, "weights", "y8", "c3"))
expect_equal(fit_ss3@fitmeasures$df, fit_ss2@fitmeasures$df - 1)
expect_true(fit_ss3@fitmeasures$chisq < fit_ss2@fitmeasures$chisq)
expect_equal(sum(fit_ss3@parameters$matrix == "weights" & fit_ss3@parameters$fixed & fit_ss3@parameters$est == 1), 6)

if (at_home()){
  # The R pipeline gives the same solution:
  fit_ss2R <- runmodel(usecpp(cvm(dat, lambda = lambda, weights = weights, beta = beta,
                                  sigma_zeta = sigma_zeta, latents = latents, sumscores = TRUE), FALSE))
  expect_true(abs(fit_ss2R@fitmeasures$chisq - fit_ss2@fitmeasures$chisq) < 1e-6)
  expect_true(max(abs(fit_ss2R@parameters$est - fit_ss2@parameters$est)) < 1e-6)

  # Comparison with lavaan (all weights fixed to 1):
  if (requireNamespace("lavaan", quietly = TRUE)){
    lav_ss <- ' f1 =~ y1 + y2 + y3
                f2 =~ y4 + y5 + y6
                c3 <~ 1*y7 + 1*y8 + 1*y9
                c4 <~ 1*y10 + 1*y11 + 1*y12 + 1*y13
                c4 ~ f1 + c3
                f2 ~ f1 + c3 + c4
                f1 ~~ c3 '
    lav_fit_ss <- tryCatch(lavaan::sem(lav_ss, data = dat[, nm], meanstructure = TRUE, likelihood = "normal",
                                       composites.cov = "free"),
                           error = function(e) NULL)
    if (!is.null(lav_fit_ss) && isTRUE(lavaan::lavInspect(lav_fit_ss, "converged"))){
      expect_equal(as.numeric(lavaan::fitMeasures(lav_fit_ss, "df")), fit_ss2@fitmeasures$df)
      expect_true(abs(as.numeric(lavaan::fitMeasures(lav_fit_ss, "chisq")) - fit_ss2@fitmeasures$chisq) < 5e-2)
      pe_ss <- lavaan::parameterEstimates(lav_fit_ss)
      b_ss <- fit_ss2@parameters[fit_ss2@parameters$matrix == "beta" & !fit_ss2@parameters$fixed, ]
      b_lav <- pe_ss[pe_ss$op == "~", ]
      b_lav <- b_lav[match(paste(b_ss$var1, b_ss$var2), paste(b_lav$lhs, b_lav$rhs)), ]
      expect_true(max(abs(b_ss$est - b_lav$est)) < 1e-2)
      expect_true(max(abs(b_ss$se - b_lav$se)) < 1e-2)
    }
  }
}

## ---- 9. Longer-running checks ----
if (at_home()){

  # Penalized ML: identical in the R pipeline, PFIML with missing data, and
  # automatic selection of the penalty:
  fit_pml_R <- suppressWarnings(suppressMessages(runmodel(usecpp(mod_pml, FALSE))))
  expect_true(max(abs(fit_pml_R@parameters$est - fit_pml@parameters$est)) < 1e-8)
  dat_mis <- dat[, nm[7:13]]
  set.seed(5); dat_mis[cbind(sample(nrow(dat_mis), 100), sample(7, 100, TRUE))] <- NA
  fit_pfiml <- suppressWarnings(suppressMessages(runmodel(pml_model(data = dat_mis, estimator = "PML", penalty_lambda = 0.3))))
  expect_equal(fit_pfiml@estimator, "PFIML")
  expect_equal(fit_pfiml@optim$convergence, 0)
  expect_equal(suppressWarnings(suppressMessages(refit(fit_pfiml)))@estimator, "FIML")
  fit_auto <- suppressWarnings(suppressMessages(runmodel(pml_model(estimator = "PML"))))
  expect_true(is.finite(fit_auto@penalty$lambda) && fit_auto@penalty$lambda > 0)
  # The model returned by the search is the penalized fit at the selected
  # penalty (the penalty must not be lost in the grid search):
  fit_fixed <- suppressWarnings(suppressMessages(runmodel(
    pml_model(estimator = "PML", penalty_lambda = fit_auto@penalty$lambda))))
  sel_edges <- fit_auto@parameters$matrix == "omega_c" & !fit_auto@parameters$fixed
  expect_true(max(abs(fit_auto@parameters$est[sel_edges] - fit_fixed@parameters$est[sel_edges])) < 1e-2)
  # At the largest penalty of the grid all penalized edges are zero:
  lam_max <- psychonetrics:::compute_lambda_max(pml_model(estimator = "PML"))
  fit_max <- suppressWarnings(suppressMessages(find_penalized_lambda(pml_model(estimator = "PML"), nLambda = 2,
                                                                     lambda_min_ratio = 0.999)))
  expect_true(all(abs(fit_max@parameters$est[sel_edges]) < 1e-3))

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
      expect_true(abs(as.numeric(lavaan::fitMeasures(lav_fit, "chisq")) - fitC@fitmeasures$chisq) < 5e-2)
      expect_true(abs(as.numeric(lavaan::fitMeasures(lav_fit, "logl")) - fitC@fitmeasures$logl) < 5e-2)
      pe <- lavaan::parameterEstimates(lav_fit)
      w_psy <- fitC@parameters[fitC@parameters$matrix == "weights" & fitC@parameters$est != 0, ]
      expect_true(max(abs(w_psy$est - pe$est[pe$op == "<~"])) < 1e-2)
      expect_true(max(abs(w_psy$se - pe$se[pe$op == "<~"]), na.rm = TRUE) < 1e-2)
      expect_true(max(abs(getmatrix(fitC, "sigma") - lavaan::lavInspect(lav_fit, "implied")$cov[nm, nm])) < 5e-2)
    }
  }
}
