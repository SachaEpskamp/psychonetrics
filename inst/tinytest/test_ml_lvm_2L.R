# Tests for the sufficient-statistics two-level ML estimator of ml_lvm()
# (estimator = "ML", distribution "TwoLevelGaussian"), added in 0.15.31.
# Phase 2 added the analytic expected Fisher information and full C++ twins
# (fit, gradient, expected Hessian, model Jacobian, prepare); the C++ path is
# now the default, and R-vs-C++ twin agreement is tested here as well.
#
# Fast deterministic checks run always (CRAN); the full lavaan comparisons
# (converged solutions, two-group, equality constraints, dyadic default) run
# under at_home() only.

suppressMessages(library(psychonetrics))
suppressMessages(library(lavaan))

# --- simulate two-level data ---
simdat_2l <- function(seed, J, sizes, p = 6, gshift = 0){
  set.seed(seed)
  nj <- if (length(sizes) == 1) rep(sizes, J) else sample(sizes, J, replace = TRUE)
  Lam <- matrix(0, p, 2); Lam[1:3,1] <- c(1,.8,.7); Lam[4:6,2] <- c(1,.9,.6)
  PsiW <- matrix(c(1,.3,.3,1),2); PsiB <- matrix(c(.5,.2,.2,.4),2)
  ThW <- diag(0.5, p); ThB <- diag(0.15, p)
  nu <- seq(0.5, 3, length = p) + gshift
  SigB <- Lam %*% PsiB %*% t(Lam) + ThB
  SigW <- Lam %*% PsiW %*% t(Lam) + ThW
  Y <- do.call(rbind, lapply(1:J, function(j){
    bj <- MASS::mvrnorm(1, rep(0,p), SigB)
    t(replicate(nj[j], nu + bj + MASS::mvrnorm(1, rep(0,p), SigW)))
  }))
  dat <- data.frame(Y); names(dat) <- paste0("y", 1:p)
  dat$cl <- rep(1:J, nj)
  dat
}

# Clean-room two-level log-likelihood from sufficient statistics (used as an
# independent oracle for the fit function):
cleanroom_logl_2l <- function(Y, cluster, mu, SW, SB){
  Y <- as.matrix(Y); N <- nrow(Y); p <- ncol(Y)
  cl <- as.integer(factor(cluster)); J <- max(cl)
  ns <- tabulate(cl)
  Ybar <- rowsum(Y, cl) / ns
  S.PW <- crossprod(Y - Ybar[cl,,drop=FALSE]) / (N - J)
  iSW <- solve(SW)
  m2ll <- (N - J) * (as.numeric(determinant(SW, TRUE)$modulus) + sum(iSW * S.PW))
  for (s in sort(unique(ns))){
    idx <- which(ns == s); m <- length(idx)
    Yb <- Ybar[idx,,drop=FALSE]
    mb <- colMeans(Yb)
    Cb <- if (m > 1) crossprod(sweep(Yb,2,mb)) / m else matrix(0,p,p)
    Sj <- SW + s * SB
    iSj <- solve(Sj)
    A <- Cb + tcrossprod(mb - mu)
    m2ll <- m2ll + m * (as.numeric(determinant(Sj, TRUE)$modulus) + s * sum(iSj * A))
  }
  -0.5 * (m2ll + N * p * log(2*pi))
}

p <- 6
lambda <- matrix(0, p, 2); lambda[1:3,1] <- 1; lambda[4:6,2] <- 1
# lavaan model with loadings constrained equal across levels (psychonetrics
# shares lambda across levels):
lmod <- '
level: 1
 fw1 =~ y1 + a2*y2 + a3*y3
 fw2 =~ y4 + a5*y5 + a6*y6
level: 2
 fb1 =~ y1 + a2*y2 + a3*y3
 fb2 =~ y4 + a5*y5 + a6*y6
'

dat <- simdat_2l(123, J = 50, sizes = 4:8)

## ---- estimator selection and model setup ----
modML <- ml_lvm(dat, lambda = lambda, clusters = "cl", estimator = "ML")
expect_equal(modML@estimator, "ML")
expect_equal(modML@distribution, "TwoLevelGaussian")
# Since Phase 2 (analytic information + C++ twins) the C++ path is the
# default, exactly as for the other model families:
expect_true(modML@cpp)
# nobs are clusters:
expect_equal(modML@sample@groups$nobs, 50L)

# default heuristic: complete data + max cluster size > 5 -> ML:
mod_def <- ml_lvm(dat, lambda = lambda, clusters = "cl")
expect_equal(mod_def@estimator, "ML")

# explicit FIML still works and is the Gaussian distribution:
modF <- ml_lvm(dat, lambda = lambda, clusters = "cl", estimator = "FIML")
expect_equal(modF@estimator, "FIML")
expect_equal(modF@distribution, "Gaussian")

# ML with within-cluster missing data is now SUPPORTED (Phase 4): it builds a
# TwoLevelGaussian model on the R path (no C++ twin for the missing-data
# likelihood) instead of erroring:
datna <- dat; datna$y1[3] <- NA
mod_na_ML <- suppressWarnings(ml_lvm(datna, lambda = lambda, clusters = "cl", estimator = "ML"))
expect_equal(mod_na_ML@estimator, "ML")
expect_equal(mod_na_ML@distribution, "TwoLevelGaussian")
expect_false(mod_na_ML@cpp)   # forced R path for missing-data ML
# the two-level statistics carry the missing-data flag:
expect_true(any(vapply(psychonetrics:::get_twolevel_stats(mod_na_ML@sample),
                       function(s) isTRUE(s$missing), logical(1))))
# default with missing data stays on FIML (conservative):
mod_na <- suppressWarnings(suppressMessages(ml_lvm(datna, lambda = lambda, clusters = "cl")))
expect_equal(mod_na@estimator, "FIML")

# default heuristic when the largest cluster has <= 5 units -> FIML (complete
# data, but clusters too small for the sufficient-statistics ML to pay off):
dat5 <- simdat_2l(11, J = 40, sizes = 5)         # all clusters size 5
expect_equal(ml_lvm(dat5, lambda = lambda, clusters = "cl")@estimator, "FIML")
dat25 <- simdat_2l(12, J = 40, sizes = 2:5)      # max size 5
expect_equal(ml_lvm(dat25, lambda = lambda, clusters = "cl")@estimator, "FIML")

# verbose = TRUE announces the chosen estimator and the experimental note;
# verbose = FALSE is silent (no experimental note). The experimental note fires
# at most once per session, so collect messages from a single verbose call:
msgs <- character(0)
withCallingHandlers(
  ml_lvm(dat, lambda = lambda, clusters = "cl", verbose = TRUE),
  message = function(m){ msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
expect_true(any(grepl("Using estimator = 'ML'", msgs)))
# verbose = FALSE: no estimator-choice or experimental messages:
msgs_quiet <- character(0)
withCallingHandlers(
  ml_lvm(dat, lambda = lambda, clusters = "cl", verbose = FALSE),
  message = function(m){ msgs_quiet <<- c(msgs_quiet, conditionMessage(m)); invokeRestart("muffleMessage") })
expect_false(any(grepl("Using estimator|experimental", msgs_quiet)))

## ---- (d) 2L objective equals the FIML objective at identical parameters ----
xs <- psychonetrics:::parVector(modML)
expect_equal(xs, psychonetrics:::parVector(modF))
fML <- psychonetrics:::psychonetrics_fitfunction(xs, modML)
fF  <- psychonetrics:::psychonetrics_fitfunction(xs, modF)
expect_true(abs(fML - fF) < 1e-10)
# and at a perturbed point:
set.seed(1); xr <- xs + runif(length(xs), -0.05, 0.05)
expect_true(abs(psychonetrics:::psychonetrics_fitfunction(xr, modML) -
                psychonetrics:::psychonetrics_fitfunction(xr, modF)) < 1e-10)

## ---- (a) log-likelihood equals the clean-room reference ----
prep <- psychonetrics:::prepareModel(xs, modML)
gm <- prep$groupModels[[1]]
ll_ref <- cleanroom_logl_2l(dat[,1:p], dat$cl, as.vector(gm$mu),
                            as.matrix(gm$sigma_within), as.matrix(gm$sigma_between))
ll_pn <- psychonetrics:::psychonetrics_logLikelihood(modML)
expect_true(abs(ll_pn - ll_ref) < 1e-10)

## ---- (c) analytic gradient equals numDeriv ----
gA <- psychonetrics:::psychonetrics_gradient(xs, modML)
gN <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, modML), xs)
expect_true(max(abs(gA - gN)) < 1e-6)

## ---- estimator cannot be switched after creation ----
expect_error(setestimator(modML, "FIML"), pattern = "cannot be switched")
expect_error(setestimator(modF, "ML"), pattern = "cannot be switched")

## ---- usecpp toggles the code path for two-level ML models ----
modML_R <- usecpp(modML, FALSE)
expect_false(modML_R@cpp)
expect_true(usecpp(modML_R, TRUE)@cpp)

## ---- Phase 2: R vs C++ twins (fit, gradient, expected Hessian, Jacobian) ----
# At the start values and at 5 random parameter vectors, on an unbalanced
# one-group model. Further configurations (balanced, two groups, at the
# optimum) are tested under at_home() below.
check_twins_2L <- function(mod, points){
  modR <- usecpp(mod, FALSE)
  modC <- usecpp(mod, TRUE)
  for (x in points){
    # Fit function:
    fR <- psychonetrics:::psychonetrics_fitfunction(x, modR)
    fC <- psychonetrics:::psychonetrics_fitfunction_cpp(x, modC)
    expect_true(abs(fR - fC) <= 1e-12 * max(1, abs(fR)))
    # Gradient:
    gR <- psychonetrics:::psychonetrics_gradient(x, modR)
    gC <- psychonetrics:::psychonetrics_gradient_cpp(x, modC)
    expect_true(max(abs(gR - gC)) <= 1e-12 * max(1, max(abs(gR))))
    # Estimator-level expected Hessian and model Jacobian on identical preps:
    prepR <- psychonetrics:::prepareModel(x, modR)
    prepC <- psychonetrics:::prepareModel(x, modC)
    hR <- as.matrix(psychonetrics:::expected_hessian_Gauss2L(prepR))
    hC <- as.matrix(psychonetrics:::expected_hessian_Gauss2L_cpp(prepC))
    expect_true(max(abs(hR - hC)) <= 1e-12 * max(1, max(abs(hR))))
    jR <- as.matrix(psychonetrics:::d_phi_theta_ml_lvm(prepR))
    jC <- as.matrix(psychonetrics:::d_phi_theta_ml_lvm_cpp(prepC))
    expect_true(max(abs(jR - jC)) <= 1e-12 * max(1, max(abs(jR))))
  }
}
xs2 <- psychonetrics:::parVector(modML)
set.seed(99)
pts <- c(list(xs2), lapply(1:5, function(i) xs2 + runif(length(xs2), -0.05, 0.05)))
check_twins_2L(modML, pts)

## ---- Phase 2: twins for non-default matrix types (augmentation paths) ----
modML_ggm <- ml_lvm(dat, lambda = lambda, clusters = "cl", estimator = "ML",
                    within_latent = "ggm", between_latent = "prec",
                    within_residual = "chol", between_residual = "ggm")
xs_g <- psychonetrics:::parVector(modML_ggm)
check_twins_2L(modML_ggm, list(xs_g))

## ---- Phase 2: analytic expected information ----
# vs the numeric (numDeriv) information of Phase 1:
Ia <- as.matrix(psychonetrics:::psychonetrics_FisherInformation(usecpp(modML, FALSE)))
In <- as.matrix(psychonetrics:::numeric_FisherInformation(usecpp(modML, FALSE)))
expect_true(max(abs(Ia - In)) / max(abs(In)) <= 1e-6)
# R vs C++ assembly of the full information matrix:
Ic <- as.matrix(psychonetrics:::psychonetrics_FisherInformation_cpp(usecpp(modML, TRUE)))
expect_true(max(abs(Ia - Ic)) <= 1e-10 * max(1, max(abs(Ia))))

## ---- Phase 2: phi-level information vs the lavaan kernel ----
# I(phi) with phi = (mu, vech Sigma_W, vech Sigma_B): psychonetrics' expected
# Hessian is twice the per-cluster unit information. lavaan's kernel orders
# blocks as (Mu.W, vech Sigma.W, Mu.B, vech Sigma.B); for ml_lvm all
# variables are both-level so the mu-gradient lives in lavaan's Mu.B block.
# lavaan >= 0.7 renamed this internal kernel and changed its interface, so the
# cross-check only runs where the 0.6 kernel is available; the analytic-vs-
# numeric information check above validates the same quantity independently:
lav_kernel <- get0("lav_mvnorm_cluster_information_expected",
                   envir = asNamespace("lavaan"), inherits = FALSE)
if (!is.null(lav_kernel)) {
  gm2 <- psychonetrics:::prepareModel(xs2, modML)$groupModels[[1]]
  fit0 <- sem(lmod, data = dat, cluster = "cl", do.fit = FALSE,
              se = "none", test = "none", baseline = FALSE)
  Lp0 <- fit0@Data@Lp[[1]]
  I_lav <- lav_kernel(
    Lp = Lp0, Mu.W = rep(0, p), Sigma.W = as.matrix(gm2$sigma_within),
    Mu.B = as.vector(gm2$mu), Sigma.B = as.matrix(gm2$sigma_between))
  k2 <- p * (p + 1) / 2
  perm <- c(p + k2 + 1:p, p + 1:k2, p + k2 + p + 1:k2)
  I_pn <- 0.5 * as.matrix(psychonetrics:::expected_hessian_Gauss2L(
    psychonetrics:::prepareModel(xs2, usecpp(modML, FALSE))))
  expect_true(max(abs(I_pn - I_lav[perm, perm])) <= 1e-10 * max(1, max(abs(I_lav))))
}

## ============================================================================
## Phase 4: within-cluster MISSING-DATA two-level ML. Fast deterministic checks
## (CRAN); full lavaan agreement runs under at_home() below.
## ============================================================================

# Impose seeded MCAR within-cluster missingness, never a fully-missing row:
impose_mcar <- function(dat, p, frac, seed){
  set.seed(seed)
  M <- matrix(runif(nrow(dat) * p) < frac, nrow(dat), p)
  for (i in which(rowSums(M) == p)) M[i, sample(p, 1)] <- FALSE
  for (k in 1:p) dat[[paste0("y", k)]][M[, k]] <- NA
  dat
}

# (m1) The missing-data likelihood reduces EXACTLY to the complete-data
# sufficient-statistics likelihood when no values are missing. Build a model on
# complete data and replace its two-level statistics by the missing-data
# structure (which on complete data has a single, fully-observed pattern), then
# compare the fit and gradient at the same parameter vector on the R path:
modR_ss <- usecpp(ml_lvm(dat, lambda = lambda, clusters = "cl", estimator = "ML"), FALSE)
ts_miss <- lapply(seq_len(nrow(modR_ss@sample@groups)), function(g){
  psychonetrics:::twolevel_missing_statistics(as.matrix(dat[, 1:p]), dat$cl)
})
modR_miss <- modR_ss
modR_miss@sample@twolevel <- ts_miss
xss <- psychonetrics:::parVector(modR_ss)
expect_true(abs(psychonetrics:::psychonetrics_fitfunction(xss, modR_ss) -
                psychonetrics:::psychonetrics_fitfunction(xss, modR_miss)) <= 1e-8)
set.seed(4); xrr <- xss + runif(length(xss), -0.03, 0.03)
expect_true(abs(psychonetrics:::psychonetrics_fitfunction(xrr, modR_ss) -
                psychonetrics:::psychonetrics_fitfunction(xrr, modR_miss)) <= 1e-8)
expect_true(max(abs(psychonetrics:::psychonetrics_gradient(xrr, modR_ss) -
                    psychonetrics:::psychonetrics_gradient(xrr, modR_miss))) <= 1e-8)

# (m2) analytic gradient of the missing-data likelihood equals numDeriv, at the
# start values and at a perturbed point (small balanced design, runs fast):
datm_s <- impose_mcar(simdat_2l(202, J = 40, sizes = 6), p, 0.15, 11)
modm_s <- suppressWarnings(ml_lvm(datm_s, lambda = lambda, clusters = "cl", estimator = "ML"))
xm <- psychonetrics:::parVector(modm_s)
gAm <- psychonetrics:::psychonetrics_gradient(xm, modm_s)
gNm <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, modm_s), xm)
expect_true(max(abs(gAm - gNm)) <= 1e-6)
set.seed(5); xm2 <- xm + runif(length(xm), -0.04, 0.04)
gAm2 <- psychonetrics:::psychonetrics_gradient(xm2, modm_s)
gNm2 <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, modm_s), xm2)
expect_true(max(abs(gAm2 - gNm2)) <= 1e-6)

# (m3) the missing-data log-likelihood (with 2*pi) matches an independent
# clean-room evaluation of the same per-pattern formula:
cleanroom_logl_2l_missing <- function(Y, cluster, mu, SW, SB){
  Y <- as.matrix(Y); p <- ncol(Y)
  cl <- as.integer(factor(cluster)); J <- max(cl)
  R <- !is.na(Y); SW.inv <- solve(SW)
  SW.ld <- as.numeric(determinant(SW, TRUE)$modulus)
  symupd <- function(na){
    if (!length(na)) return(list(inv = SW.inv, ld = SW.ld))
    H <- SW.inv[na, na, drop = FALSE]; A <- SW.inv[na, -na, drop = FALSE]
    list(inv = SW.inv[-na, -na, drop = FALSE] - crossprod(A, solve(H, A)),
         ld = SW.ld + as.numeric(determinant(H, TRUE)$modulus))
  }
  Yc <- sweep(Y, 2, mu); PIJ <- matrix(0, nrow(Y), p)
  Wlog <- 0; Alist <- rep(list(matrix(0, p, p)), J)
  patkey <- apply(R, 1, function(r) paste(which(r), collapse = ","))
  for (k in unique(patkey)){
    rows <- which(patkey == k); o <- which(R[rows[1], ]); na <- which(!R[rows[1], ])
    u <- symupd(na); Wlog <- Wlog + u$ld * length(rows)
    PIJ[rows, o] <- Yc[rows, o, drop = FALSE] %*% u$inv
    Wf <- matrix(0, p, p); Wf[o, o] <- u$inv
    for (j in unique(cl[rows])) Alist[[j]] <- Alist[[j]] + Wf * sum(cl[rows] == j)
  }
  qa <- sum(PIJ * Yc, na.rm = TRUE)
  PJ <- rowsum(PIJ, cl, reorder = TRUE); qb <- 0; ld <- 0
  for (j in 1:J){
    M <- SB %*% Alist[[j]]; diag(M) <- diag(M) + 1
    ld <- ld + as.numeric(determinant(M, TRUE)$modulus)
    qb <- qb + sum(PJ[j, ] * solve(M, drop(SB %*% PJ[j, ])))
  }
  -0.5 * ((qa - qb) + (Wlog + ld) + sum(R) * log(2 * pi))
}
gmm <- psychonetrics:::prepareModel(xm, modm_s)$groupModels[[1]]
ll_ref_m <- cleanroom_logl_2l_missing(datm_s[, 1:p], datm_s$cl, as.vector(gmm$mu),
                                      as.matrix(gmm$sigma_within), as.matrix(gmm$sigma_between))
ll_pn_m <- psychonetrics:::psychonetrics_logLikelihood(modm_s)
expect_true(abs(ll_pn_m - ll_ref_m) < 1e-8)

## ---- full lavaan comparisons (slow) ----
if (at_home()){

  ## (b) converged solution vs lavaan, unbalanced:
  datu <- simdat_2l(321, J = 100, sizes = 2:25)
  modu <- ml_lvm(datu, lambda = lambda, clusters = "cl", estimator = "ML")
  modu@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
  modu <- runmodel(modu, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE)
  fitu <- sem(lmod, data = datu, cluster = "cl", se = "none", test = "none", baseline = FALSE)
  expect_true(abs(modu@fitmeasures$logl - as.numeric(logLik(fitu))) < 1e-5)
  # estimates via implied moments:
  xo <- psychonetrics:::parVector(modu)
  gmo <- psychonetrics:::prepareModel(xo, modu)$groupModels[[1]]
  impu <- lavInspect(fitu, "implied")
  expect_true(max(abs(as.matrix(gmo$sigma_within) - as.matrix(impu$within$cov))) < 1e-4)
  expect_true(max(abs(as.matrix(gmo$sigma_between) - as.matrix(impu$cl$cov))) < 1e-4)
  expect_true(max(abs(as.vector(gmo$mu) -
                      (as.numeric(impu$within$mean) + as.numeric(impu$cl$mean)))) < 1e-4)
  # gradient at the optimum:
  gAo <- psychonetrics:::psychonetrics_gradient(xo, modu)
  gNo <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, modu), xo)
  expect_true(max(abs(gAo - gNo)) < 1e-6)

  ## R vs C++ twins at the (unbalanced) optimum and on a balanced model:
  modu_opt <- psychonetrics:::updateModel(xo, modu)
  set.seed(77)
  pts_opt <- c(list(xo), lapply(1:5, function(i) xo + runif(length(xo), -0.05, 0.05)))
  check_twins_2L(modu_opt, pts_opt)

  ## (b) balanced, including df/chisq and SEs:
  datb <- simdat_2l(42, J = 100, sizes = 10)
  modb <- ml_lvm(datb, lambda = lambda, clusters = "cl", estimator = "ML")
  modb@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
  xb <- psychonetrics:::parVector(modb)
  set.seed(78)
  check_twins_2L(modb, c(list(xb), lapply(1:5, function(i) xb + runif(length(xb), -0.05, 0.05))))
  modb <- runmodel(modb, addMIs = FALSE)
  fitb <- sem(lmod, data = datb, cluster = "cl")
  expect_true(abs(modb@fitmeasures$logl - as.numeric(logLik(fitb))) < 1e-5)
  expect_equal(as.numeric(modb@fitmeasures$df), as.numeric(fitMeasures(fitb, "df")))
  expect_true(abs(modb@fitmeasures$chisq - fitMeasures(fitb, "chisq")) < 1e-2)
  # SEs (analytic expected Fisher information):
  prm <- modb@parameters
  expect_false(any(is.na(prm$se[!prm$fixed])))

  # SEs match lavaan with information = "expected" (lavaan's multilevel
  # DEFAULT is the observed information, which differs by a few percent):
  fitb_exp <- sem(lmod, data = datb, cluster = "cl", information = "expected")
  pe <- parameterEstimates(fitb_exp)
  freeb <- prm[!prm$fixed & !duplicated(prm$par), ]
  match_se <- sapply(seq_len(nrow(freeb)), function(i){
    d <- abs(pe$est - freeb$est[i])
    j <- which.min(d)
    if (d[j] < 1e-3) pe$se[j] else NA
  })
  okb <- !is.na(match_se) & match_se > 0
  expect_true(sum(okb) >= nrow(freeb) - 2) # nearly all matched
  expect_true(max(abs(freeb$se[okb] - match_se[okb]) / match_se[okb]) < 1e-4)

  ## full runmodel: C++ path identical to the R path:
  modbR <- usecpp(ml_lvm(datb, lambda = lambda, clusters = "cl", estimator = "ML"), FALSE)
  modbR@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                          iter.max = 1000, eval.max = 2000))
  modbR <- runmodel(modbR, addMIs = FALSE)
  expect_true(abs(modb@fitmeasures$logl - modbR@fitmeasures$logl) < 1e-8)
  expect_true(max(abs(modb@parameters$est - modbR@parameters$est)) < 1e-8)
  expect_true(max(abs(modb@parameters$se - modbR@parameters$se), na.rm = TRUE) < 1e-8)

  ## (e) two-group model. lavaan 0.6-21 cannot fit cluster= and group=
  ## together, but the unconstrained multigroup model equals separate fits:
  d1 <- simdat_2l(7, 60, 5:12); d1$gr <- "g1"
  d2 <- simdat_2l(8, 70, 5:12, gshift = 0.3); d2$gr <- "g2"; d2$cl <- d2$cl + 1000
  datg <- rbind(d1, d2)
  modg <- ml_lvm(datg, lambda = lambda, clusters = "cl", groups = "gr", estimator = "ML")
  modg@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
  # R vs C++ twins on the two-group model:
  xg <- psychonetrics:::parVector(modg)
  set.seed(79)
  check_twins_2L(modg, c(list(xg), lapply(1:5, function(i) xg + runif(length(xg), -0.05, 0.05))))
  modg <- runmodel(modg, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE)
  fit1 <- sem(lmod, data = d1, cluster = "cl", se = "none", test = "none", baseline = FALSE)
  fit2 <- sem(lmod, data = d2, cluster = "cl", se = "none", test = "none", baseline = FALSE)
  expect_true(abs(modg@fitmeasures$logl -
                  (as.numeric(logLik(fit1)) + as.numeric(logLik(fit2)))) < 1e-5)

  ## (f) equality-constrained two-group model:
  modge <- ml_lvm(datg, lambda = lambda, clusters = "cl", groups = "gr",
                  estimator = "ML", equal = "lambda")
  modge <- runmodel(modge, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE)
  prge <- modge@parameters
  expect_identical(prge$est[prge$matrix == "lambda" & prge$group_id == 1],
                   prge$est[prge$matrix == "lambda" & prge$group_id == 2])
  expect_true(modge@fitmeasures$df > modg@fitmeasures$df)
  expect_true(modge@fitmeasures$logl <= modg@fitmeasures$logl + 1e-8)

  ## (h) dyadic data: default heuristic picks FIML and the model runs:
  datd <- simdat_2l(9, J = 100, sizes = 2)
  modd <- suppressMessages(ml_lvm(datd, lambda = lambda, clusters = "cl"))
  expect_equal(modd@estimator, "FIML")
  modd <- suppressWarnings(runmodel(modd, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE))
  expect_true(modd@computed)

  ## (1) saturated/baseline reference models and fit measures vs lavaan.
  ## The saturated and baseline models are themselves two-level models that
  ## runmodel() optimizes numerically; these must reach the same maximum as
  ## lavaan's EM-based h1/baseline. Use the balanced model fitted above (modb,
  ## fitb) so the comparison is on a converged, well-specified model.
  lavb <- fitMeasures(fitb)
  # Saturated (unrestricted) log-likelihood: psychonetrics optimizes it, lavaan
  # uses EM (default tol ~1e-4 in older lavaan), so allow 1e-3:
  expect_true(abs(modb@fitmeasures$unrestricted.logl -
                  as.numeric(lavb["unrestricted.logl"])) < 1e-3)
  # Baseline (independence) log-likelihood vs lavaan's independence model.
  # lavaan's fitMeasures() does not always expose baseline.logl for two-level
  # models, so refit the independence model explicitly:
  basll_lav <- as.numeric(logLik(lavaan:::lav_object_independence(fitb)))
  expect_true(abs(modb@fitmeasures$baseline.logl - basll_lav) < 1e-3)
  # Model chi-square and degrees of freedom: chisq close, df EXACTLY equal:
  expect_true(abs(modb@fitmeasures$chisq - as.numeric(lavb["chisq"])) < 1e-3)
  expect_equal(as.numeric(modb@fitmeasures$df), as.numeric(lavb["df"]))
  # Baseline chi-square and df (drives CFI/TLI): df EXACTLY equal:
  expect_true(abs(modb@fitmeasures$baseline.chisq - as.numeric(lavb["baseline.chisq"])) < 1e-2)
  expect_equal(as.numeric(modb@fitmeasures$baseline.df), as.numeric(lavb["baseline.df"]))
  # CFI/TLI finite, plausible, and agree with lavaan (same baseline convention):
  expect_true(is.finite(modb@fitmeasures$cfi) && modb@fitmeasures$cfi > 0.9 &&
              modb@fitmeasures$cfi <= 1 + 1e-8)
  expect_true(is.finite(modb@fitmeasures$tli))
  expect_true(abs(modb@fitmeasures$cfi - as.numeric(lavb["cfi"])) < 1e-3)
  expect_true(abs(modb@fitmeasures$tli - as.numeric(lavb["tli"])) < 1e-3)

  ## (1b) RMSEA n-convention: psychonetrics uses n = #clusters (J), lavaan uses
  ## n = #units (N). On a MISSPECIFIED model (nonzero RMSEA) with chi-square and
  ## df otherwise agreeing, psychonetrics' RMSEA is larger by EXACTLY sqrt(N/J).
  ## Fit a 1-factor model (the 2-factor data misfit it) with loadings shared
  ## across levels, and a matching lavaan model with equal level-1/level-2
  ## loadings:
  lambda1 <- matrix(1, p, 1)
  lmod1 <- '
  level: 1
   fw =~ y1 + L2*y2 + L3*y3 + L4*y4 + L5*y5 + L6*y6
  level: 2
   fb =~ y1 + L2*y2 + L3*y3 + L4*y4 + L5*y5 + L6*y6
  '
  modm <- ml_lvm(datb, lambda = lambda1, clusters = "cl", estimator = "ML")
  modm@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
  modm <- runmodel(modm, addMIs = FALSE)
  fitm <- sem(lmod1, data = datb, cluster = "cl")
  Nb <- nrow(datb); Jb <- length(unique(datb$cl))
  # chi-square and df agree with lavaan (well-matched parameterization):
  expect_true(abs(modm@fitmeasures$chisq - fitMeasures(fitm, "chisq")) < 1e-2)
  expect_equal(as.numeric(modm@fitmeasures$df), as.numeric(fitMeasures(fitm, "df")))
  # nonzero RMSEA, and the ratio to lavaan's is exactly sqrt(N/J):
  expect_true(modm@fitmeasures$rmsea > 0.05)
  expect_true(abs(modm@fitmeasures$rmsea /
                  as.numeric(fitMeasures(fitm, "rmsea")) - sqrt(Nb/Jb)) < 1e-4)
  # internal consistency: ML's misspecified-model logl equals FIML's (the two
  # estimators optimize the same objective), even though both differ from
  # lavaan when the model is misspecified and clusters are unbalanced:
  modmF <- ml_lvm(datb, lambda = lambda1, clusters = "cl", estimator = "FIML")
  modmF@optim.args <- modm@optim.args
  modmF <- runmodel(modmF, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE)
  expect_true(abs(modm@fitmeasures$logl - modmF@fitmeasures$logl) < 1e-4)

  ## (3) performance: a full runmodel under "ML" on J = 100 clusters of size
  ## 18-25 with p = 6 is fast. The sufficient-statistics cost does not grow with
  ## cluster size, so this is well under a second on CI; use a generous bound:
  datp <- simdat_2l(2024, J = 100, sizes = 18:25)
  tp <- system.time(
    modp <- runmodel(ml_lvm(datp, lambda = lambda, clusters = "cl", estimator = "ML"))
  )[["elapsed"]]
  expect_true(modp@computed)
  expect_true(is.finite(modp@fitmeasures$cfi))
  expect_true(tp < 5)

  ## ==========================================================================
  ## Phase 4 (slow): within-cluster MISSING-DATA two-level ML vs lavaan with
  ## missing = "ml", for an unbalanced and a balanced design.
  ## ==========================================================================

  ## (m4) UNBALANCED, 15% MCAR: logl, implied moments and gradient vs lavaan:
  datmu <- impose_mcar(simdat_2l(321, J = 80, sizes = 5:15), p, 0.15, 9001)
  modmu <- suppressWarnings(ml_lvm(datmu, lambda = lambda, clusters = "cl", estimator = "ML"))
  modmu@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                          iter.max = 2000, eval.max = 3000))
  modmu <- runmodel(modmu, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE)
  fitmu <- sem(lmod, data = datmu, cluster = "cl", missing = "ml",
               se = "none", test = "none", baseline = FALSE)
  expect_true(abs(modmu@fitmeasures$logl - as.numeric(logLik(fitmu))) < 1e-4)
  xou <- psychonetrics:::parVector(modmu)
  gmou <- psychonetrics:::prepareModel(xou, modmu)$groupModels[[1]]
  implu <- lavInspect(fitmu, "implied")
  expect_true(max(abs(as.matrix(gmou$sigma_within) - as.matrix(implu$within$cov))) < 1e-3)
  expect_true(max(abs(as.matrix(gmou$sigma_between) - as.matrix(implu$cl$cov))) < 1e-3)
  expect_true(max(abs(as.vector(gmou$mu) -
                      (as.numeric(implu$within$mean) + as.numeric(implu$cl$mean)))) < 1e-3)
  gAou <- psychonetrics:::psychonetrics_gradient(xou, modmu)
  gNou <- numDeriv::grad(function(par) psychonetrics:::psychonetrics_fitfunction(par, modmu), xou)
  expect_true(max(abs(gAou - gNou)) < 1e-6)

  ## (m5) BALANCED, 12% MCAR: logl, df, chisq and SEs vs lavaan (observed
  ## information, which is lavaan's multilevel default and the convention used
  ## by psychonetrics' numeric Fisher information for the missing-data path):
  datmb <- impose_mcar(simdat_2l(42, J = 100, sizes = 10), p, 0.12, 9002)
  modmb <- suppressWarnings(ml_lvm(datmb, lambda = lambda, clusters = "cl", estimator = "ML"))
  modmb@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                          iter.max = 2000, eval.max = 3000))
  modmb <- runmodel(modmb, addMIs = FALSE)
  fitmb <- sem(lmod, data = datmb, cluster = "cl", missing = "ml",
               information = "observed")
  expect_true(abs(modmb@fitmeasures$logl - as.numeric(logLik(fitmb))) < 1e-4)
  expect_equal(as.numeric(modmb@fitmeasures$df), as.numeric(fitMeasures(fitmb, "df")))
  expect_true(abs(modmb@fitmeasures$chisq - fitMeasures(fitmb, "chisq")) < 1e-1)
  # SEs are finite for all free parameters and match lavaan's observed-info SEs:
  prmb <- modmb@parameters
  expect_false(any(is.na(prmb$se[!prmb$fixed])))
  peb <- parameterEstimates(fitmb)
  freebm <- prmb[!prmb$fixed & !duplicated(prmb$par), ]
  match_se_m <- sapply(seq_len(nrow(freebm)), function(i){
    d <- abs(peb$est - freebm$est[i]); j <- which.min(d)
    if (d[j] < 1e-2) peb$se[j] else NA
  })
  okm <- !is.na(match_se_m) & match_se_m > 0
  expect_true(sum(okm) >= nrow(freebm) - 2)
  expect_true(max(abs(freebm$se[okm] - match_se_m[okm]) / match_se_m[okm]) < 1e-2)
  # (Note: with missing data the two-level ML and the wide-format FIML do NOT
  # optimize the same objective -- FIML's missingness patterns are over the
  # wide-format cluster rows, not the per-unit two-level patterns -- so they
  # need not coincide, unlike the complete-data case. The relevant external
  # reference is lavaan's two-level missing-data ML, checked in (m4)/(m5).)
}
