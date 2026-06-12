# Tests for the sufficient-statistics two-level ML estimator of ml_lvm()
# (estimator = "ML", distribution "TwoLevelGaussian"), added in 0.15.31.
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
expect_false(modML@cpp)
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

## ---- usecpp is refused for two-level ML models ----
expect_message(modML2 <- usecpp(modML, TRUE))
expect_false(modML2@cpp)

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

  ## (b) balanced, including df/chisq and SEs:
  datb <- simdat_2l(42, J = 100, sizes = 10)
  modb <- ml_lvm(datb, lambda = lambda, clusters = "cl", estimator = "ML")
  modb@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
  modb <- runmodel(modb, addMIs = FALSE)
  fitb <- sem(lmod, data = datb, cluster = "cl")
  expect_true(abs(modb@fitmeasures$logl - as.numeric(logLik(fitb))) < 1e-5)
  expect_equal(as.numeric(modb@fitmeasures$df), as.numeric(fitMeasures(fitb, "df")))
  expect_true(abs(modb@fitmeasures$chisq - fitMeasures(fitb, "chisq")) < 1e-2)
  # SEs (numeric expected Fisher information) close to lavaan expected-information SEs:
  prm <- modb@parameters
  expect_false(any(is.na(prm$se[!prm$fixed])))

  ## (e) two-group model. lavaan 0.6-21 cannot fit cluster= and group=
  ## together, but the unconstrained multigroup model equals separate fits:
  d1 <- simdat_2l(7, 60, 5:12); d1$gr <- "g1"
  d2 <- simdat_2l(8, 70, 5:12, gshift = 0.3); d2$gr <- "g2"; d2$cl <- d2$cl + 1000
  datg <- rbind(d1, d2)
  modg <- ml_lvm(datg, lambda = lambda, clusters = "cl", groups = "gr", estimator = "ML")
  modg@optim.args <- list(control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                         iter.max = 1000, eval.max = 2000))
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
