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

# ML with missing data errors informatively:
datna <- dat; datna$y1[3] <- NA
expect_error(suppressWarnings(ml_lvm(datna, lambda = lambda, clusters = "cl", estimator = "ML")),
             pattern = "missing data")
# default with missing data picks FIML:
mod_na <- suppressWarnings(suppressMessages(ml_lvm(datna, lambda = lambda, clusters = "cl")))
expect_equal(mod_na@estimator, "FIML")

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
# variables are both-level so the mu-gradient lives in lavaan's Mu.B block:
gm2 <- psychonetrics:::prepareModel(xs2, modML)$groupModels[[1]]
fit0 <- sem(lmod, data = dat, cluster = "cl", do.fit = FALSE,
            se = "none", test = "none", baseline = FALSE)
Lp0 <- fit0@Data@Lp[[1]]
I_lav <- lavaan:::lav_mvnorm_cluster_information_expected(
  Lp = Lp0, Mu.W = rep(0, p), Sigma.W = as.matrix(gm2$sigma_within),
  Mu.B = as.vector(gm2$mu), Sigma.B = as.matrix(gm2$sigma_between))
k2 <- p * (p + 1) / 2
perm <- c(p + k2 + 1:p, p + 1:k2, p + k2 + p + 1:k2)
I_pn <- 0.5 * as.matrix(psychonetrics:::expected_hessian_Gauss2L(
  psychonetrics:::prepareModel(xs2, usecpp(modML, FALSE))))
expect_true(max(abs(I_pn - I_lav[perm, perm])) <= 1e-10 * max(1, max(abs(I_lav))))

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
}
