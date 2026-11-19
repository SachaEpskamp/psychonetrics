# Tests for the lavaan interoperability functions fromlavaan() and tolavaan().
#
# The fast subset (always run) exercises the conversion machinery and all of
# the guard conditions on small models. The numerical agreement with lavaan
# 0.6.21 (reference logLik/chisq/df/npar and parameter/SE differences) is
# checked under at_home() only, since it depends on lavaan being installed and
# on a full runmodel()/lavaan() fit.
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build) against lavaan 0.6.21.

suppressMessages(library(psychonetrics))
if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")
suppressMessages(library(lavaan))

data(HolzingerSwineford1939, package = "lavaan")
HS <- lavaan::HolzingerSwineford1939

Lambda9 <- matrix(0, 9, 3)
Lambda9[1:3, 1] <- Lambda9[4:6, 2] <- Lambda9[7:9, 3] <- 1

## ============================ fromlavaan() ============================
## Fast: single-group CFA, marker identification + meanstructure.
f1 <- suppressWarnings(cfa(' visual  =~ x1 + x2 + x3
                            textual =~ x4 + x5 + x6
                            speed   =~ x7 + x8 + x9 ',
                          data = HS, meanstructure = TRUE))
m1 <- fromlavaan(f1)
expect_inherits(m1, "psychonetrics")
expect_equal(m1@model, "lvm")
expect_equal(m1@types$latent, "cov")
expect_equal(m1@types$residual, "cov")
# Conversion must not fit unless run = TRUE:
expect_false(m1@computed)
# Skeleton must contain a free lambda for every indicator-loading and the three
# fixed marker loadings (value 1, identified):
lam <- m1@parameters[m1@parameters$matrix == "lambda", ]
expect_equal(sum(lam$est == 1 & lam$fixed), 3L)         # markers
expect_equal(sum(!lam$fixed), 6L)                       # freely estimated loadings

## Fast: syntax + data entry point yields an identical skeleton.
m1b <- fromlavaan(' visual  =~ x1 + x2 + x3
                    textual =~ x4 + x5 + x6
                    speed   =~ x7 + x8 + x9 ', data = HS)
expect_equal(nrow(m1b@parameters), nrow(m1@parameters))

## Fast: syntax input without data is an error.
expect_error(fromlavaan(' visual =~ x1 + x2 + x3 '), "data")

## Fast: unsupported lavaan features raise informative errors.
expect_error(fromlavaan(suppressWarnings(cfa(
  ' visual =~ x1 + x2 + x3
    d := b1 * b1
    b1 := 1 ', data = HS))), "Defined parameters")
# Categorical endogenous variable:
HSc <- HS; HSc$x1o <- ordered(cut(HS$x1, 3))
expect_error(suppressWarnings(fromlavaan(cfa(
  ' visual =~ x1o + x2 + x3 ', data = HSc, ordered = "x1o"))),
  "[Cc]ategorical|ordinal")

## ============================= tolavaan() =============================
## Fast: single-group lvm(cov/cov) -> parameter table (fit = FALSE).
mcfa <- suppressMessages(runmodel(lvm(HS, lambda = Lambda9,
  latent = "cov", residual = "cov", vars = paste0("x", 1:9),
  latents = c("visual", "textual", "speed"), identification = "loadings")))
PT <- tolavaan(mcfa, fit = FALSE)
expect_inherits(PT, "data.frame")
expect_true(all(c("lhs", "op", "rhs", "block", "group", "free",
                  "ustart", "label", "plabel", "id") %in% names(PT)))
# Loadings, (co)variances, but no mean rows (meanstructure off here by default
# only if nu absent) -- check the operator set is a subset of lavaan ops:
expect_true(all(PT$op %in% c("=~", "~~", "~", "~1", "==")))
# Free indices are a contiguous 1..k sequence on the non-fixed rows:
fr <- sort(PT$free[PT$free > 0])
expect_equal(fr, seq_len(length(fr)))

## Fast: guards fire (informative errors).
expect_error(tolavaan(42), "psychonetrics")
vc <- suppressMessages(varcov(HS, vars = paste0("x", 1:9)))
expect_error(tolavaan(vc), "lvm-family")
mggm <- suppressMessages(lvm(HS, lambda = Lambda9, latent = "ggm",
  residual = "cov", vars = paste0("x", 1:9),
  latents = c("visual", "textual", "speed"), identification = "variance"))
expect_error(tolavaan(mggm), "cov")

## ====================== at_home: lavaan agreement =====================
if (at_home()) {
  data(PoliticalDemocracy, package = "lavaan")
  PD <- lavaan::PoliticalDemocracy

  ## Map free lavaan parameters to psychonetrics rows; return est/SE diffs.
  parcmp <- function(fit, mod){
    nG <- lavInspect(fit, "ngroups")
    FREE <- lavInspect(fit, "free"); if (nG == 1) FREE <- list(FREE)
    mm <- c(lambda = "lambda", theta = "sigma_epsilon", psi = "sigma_zeta",
            beta = "beta", nu = "nu", alpha = "nu_eta")
    pt <- parTable(fit); pt <- pt[pt$free > 0, ]
    P <- mod@parameters
    out <- do.call(rbind, lapply(seq_len(nrow(pt)), function(r){
      f <- pt$free[r]
      for (g in 1:nG) for (m in names(mm)){
        Fm <- FREE[[g]][[m]]; if (is.null(Fm)) next
        w <- which(Fm == f, arr.ind = TRUE)
        if (nrow(w) > 0){
          i <- w[1, 1]; j <- if (ncol(Fm) > 1 || m %in% c("nu", "alpha")) w[1, 2] else 1
          if (m %in% c("theta", "psi") && j > i){ t <- i; i <- j; j <- t }
          idx <- which(P$matrix == mm[[m]] & P$row == i & P$col == j & P$group_id == g)
          return(data.frame(ed = abs(pt$est[r] - P$est[idx]),
                            sd = abs(pt$se[r]  - P$se[idx])))
        }
      }
      NULL
    }))
    c(est = max(out$ed), se = max(out$sd))
  }

  ## (1) HS 3-factor CFA, marker + meanstructure.
  M1 <- suppressMessages(runmodel(m1))
  expect_true(abs(M1@fitmeasures$logl - (-3737.744927)) < 1e-5)
  expect_true(abs(M1@fitmeasures$chisq - 85.305522) < 1e-4)
  expect_equal(M1@fitmeasures$df, 24)
  expect_equal(max(M1@parameters$par), 30L)
  d1 <- parcmp(f1, M1)
  expect_true(d1["est"] < 1e-5)
  expect_true(d1["se"]  < 1e-4)

  ## (2) PoliticalDemocracy SEM (beta + residual covariances + meanstructure).
  f2 <- sem('
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8 ', data = PD, meanstructure = TRUE)
  M2 <- suppressMessages(runmodel(fromlavaan(f2)))
  expect_true(abs(M2@fitmeasures$logl - (-1547.790943)) < 1e-5)
  expect_true(abs(M2@fitmeasures$chisq - 38.125218) < 1e-4)
  expect_equal(M2@fitmeasures$df, 35)
  expect_equal(max(M2@parameters$par), 42L)
  d2 <- parcmp(f2, M2)
  expect_true(d2["est"] < 1e-7)
  expect_true(d2["se"]  < 1e-7)

  ## (3) Two-group CFA with group.equal = c(loadings, intercepts).
  f3 <- suppressWarnings(cfa(' visual  =~ x1 + x2 + x3
                              textual =~ x4 + x5 + x6
                              speed   =~ x7 + x8 + x9 ',
                            data = HS, group = "school",
                            group.equal = c("loadings", "intercepts"),
                            meanstructure = TRUE))
  M3 <- suppressMessages(runmodel(fromlavaan(f3)))
  expect_true(abs(M3@fitmeasures$logl - (-3706.323256)) < 1e-4)
  expect_true(abs(M3@fitmeasures$chisq - 164.102831) < 1e-3)
  expect_equal(M3@fitmeasures$df, 60)
  expect_equal(max(M3@parameters$par), 48L)

  ## (4) std.lv identification.
  f4 <- suppressWarnings(cfa(' visual =~ x1 + x2 + x3 ',
                            data = HS, std.lv = TRUE, meanstructure = TRUE))
  M4 <- suppressMessages(runmodel(fromlavaan(f4)))
  expect_true(abs(M4@fitmeasures$logl - (-1356.977317)) < 1e-4)

  ## (5) FIML with injected missingness: logLik matches lavaan.
  set.seed(1); HSm <- HS
  for (v in paste0("x", 1:9)) HSm[sample(nrow(HSm), 20), v] <- NA
  f5 <- suppressWarnings(cfa(' visual  =~ x1 + x2 + x3
                              textual =~ x4 + x5 + x6
                              speed   =~ x7 + x8 + x9 ',
                            data = HSm, meanstructure = TRUE, missing = "fiml"))
  M5 <- suppressMessages(runmodel(fromlavaan(f5)))
  expect_true(abs(M5@fitmeasures$logl - as.numeric(logLik(f5))) < 1e-6)

  ## tolavaan(): 2-group lvm(cov/cov) with equal loadings -> lavaan.
  modA <- suppressMessages(runmodel(groupequal(lvm(HS, lambda = Lambda9,
    latent = "cov", residual = "cov", vars = paste0("x", 1:9),
    latents = c("visual", "textual", "speed"),
    groupvar = "school", identification = "loadings"), "lambda")))
  lavA <- tolavaan(modA)
  expect_equal(modA@fitmeasures$df, 54)
  expect_equal(as.integer(fitMeasures(lavA, "df")), 54L)
  # Warm start: lavaan polishes a few steps, so allow optimizer tolerance.
  expect_true(abs(modA@fitmeasures$logl - as.numeric(logLik(lavA))) < 1e-4)
  expect_true(abs(modA@fitmeasures$chisq - fitMeasures(lavA, "chisq")) < 1e-3)
  # Parameter / SE agreement:
  P <- modA@parameters
  P <- P[!P$fixed & P$matrix %in% c("lambda", "sigma_zeta", "sigma_epsilon",
                                    "beta", "nu", "nu_eta"), ]
  pt <- parTable(lavA)
  vars <- paste0("x", 1:9); lat <- c("visual", "textual", "speed")
  getlav <- function(m, r, c, g){
    trip <- switch(m,
      lambda = c(lat[c], "=~", vars[r]), sigma_zeta = c(lat[c], "~~", lat[r]),
      sigma_epsilon = c(vars[c], "~~", vars[r]), beta = c(lat[r], "~", lat[c]),
      nu = c(vars[r], "~1", ""), nu_eta = c(lat[r], "~1", ""))
    i <- which(pt$lhs == trip[1] & pt$op == trip[2] & pt$rhs == trip[3] & pt$group == g)
    c(pt$est[i], pt$se[i])
  }
  cmp <- t(mapply(getlav, P$matrix, P$row, P$col, P$group_id))
  expect_true(max(abs(P$est - cmp[, 1])) < 1e-3)
  expect_true(max(abs(P$se  - cmp[, 2])) < 1e-4)

  ## tolavaan(): fit = FALSE for this model has exactly 6 "==" rows.
  PTa <- tolavaan(modA, fit = FALSE)
  expect_equal(sum(PTa$op == "=="), 6L)

  ## tolavaan(): FIML with storedata = TRUE + missing data.
  modF <- suppressMessages(runmodel(lvm(HSm, lambda = Lambda9,
    latent = "cov", residual = "cov", vars = paste0("x", 1:9),
    latents = c("visual", "textual", "speed"), identification = "loadings",
    estimator = "FIML", storedata = TRUE)))
  lavF <- tolavaan(modF)
  # At the warm-start point lavaan's logLik matches exactly; the refit can
  # polish a few optimizer steps, so use the same tolerance as the ML case.
  expect_true(abs(modF@fitmeasures$logl - as.numeric(logLik(lavF))) < 1e-4)
  # do.fit = FALSE evaluates lavaan AT the psychonetrics estimates -> exact:
  lavF0 <- tolavaan(modF, do.fit = FALSE)
  expect_true(abs(modF@fitmeasures$logl - as.numeric(logLik(lavF0))) < 1e-6)

  ## tolavaan(): FIML without stored data is an error.
  modFns <- suppressMessages(runmodel(lvm(HSm, lambda = Lambda9,
    latent = "cov", residual = "cov", vars = paste0("x", 1:9),
    latents = c("visual", "textual", "speed"), identification = "loadings",
    estimator = "FIML", storedata = FALSE)))
  expect_error(tolavaan(modFns), "storedata")

  ## Round-trip: tolavaan(fromlavaan(fit)) reproduces the original logLik.
  lavRT <- tolavaan(suppressMessages(runmodel(fromlavaan(f1))))
  expect_true(abs(as.numeric(logLik(f1)) - as.numeric(logLik(lavRT))) < 1e-4)
}

## =========================================================================
## Regression (0.16.5): fromlavaan() on a DEFAULT (meanstructure = FALSE)
## single-group fit must not collide nu with the covariance-structure
## parameters. Previously chi-square was Inf and df was wrong.
## =========================================================================
if (at_home()){
  hsmod <- "visual =~ x1+x2+x3\ntextual =~ x4+x5+x6\nspeed =~ x7+x8+x9"

  # (a) default cfa() (no mean structure):
  f_nom <- lavaan::cfa(hsmod, HolzingerSwineford1939)
  p_nom <- suppressMessages(fromlavaan(f_nom, run = TRUE))
  expect_equal(p_nom@fitmeasures$df, as.integer(lavaan::fitMeasures(f_nom)["df"]))
  expect_true(abs(p_nom@fitmeasures$chisq -
                    as.numeric(lavaan::fitMeasures(f_nom)["chisq"])) < 1e-3)
  # no par collision: nu indices disjoint from lambda indices
  pt <- p_nom@parameters
  expect_true(length(intersect(pt$par[pt$matrix == "nu" & !pt$fixed],
                               pt$par[pt$matrix == "lambda" & !pt$fixed])) == 0)

  # (b) std.lv default:
  f_std <- lavaan::cfa("f =~ x1+x2+x3+x4", HolzingerSwineford1939, std.lv = TRUE)
  p_std <- suppressMessages(fromlavaan(f_std, run = TRUE))
  expect_true(abs(p_std@fitmeasures$chisq -
                    as.numeric(lavaan::fitMeasures(f_std)["chisq"])) < 1e-3)

  # (c) default sem() with a regression:
  f_sem <- lavaan::sem("ind60 =~ x1+x2+x3\ndem60 =~ y1+y2+y3+y4\ndem60 ~ ind60",
                       lavaan::PoliticalDemocracy)
  p_sem <- suppressMessages(fromlavaan(f_sem, run = TRUE))
  expect_equal(p_sem@fitmeasures$df, as.integer(lavaan::fitMeasures(f_sem)["df"]))
  expect_true(abs(p_sem@fitmeasures$chisq -
                    as.numeric(lavaan::fitMeasures(f_sem)["chisq"])) < 1e-3)
}


## =========================================================================
## Composite variable models: cvm() <-> lavaan composites ("<~")
## =========================================================================
vars9 <- paste0("x", 1:9)
Wc <- matrix(0, 9, 3); Wc[1:3, 1] <- 1; Wc[4:6, 2] <- 1
Lc <- matrix(0, 9, 3); Lc[7:9, 3] <- 1
Bc <- matrix(0, 3, 3); Bc[3, 1:2] <- 1
Szc <- diag(3); Szc[1, 2] <- Szc[2, 1] <- 1
cvm_mod <- cvm(HS, vars = vars9, lambda = Lc, weights = Wc, beta = Bc, sigma_zeta = Szc,
               latents = c("c1", "c2", "f"))

## Fast: the lavaan parameter table of a cvm (no lavaan fit involved).
if (utils::packageVersion("lavaan") >= "0.7.0"){
  PTc <- tolavaan(cvm_mod, fit = FALSE)
  expect_equal(sum(PTc$op == "<~"), 6L)                     # weights
  expect_equal(sum(PTc$op == "=~"), 3L)                     # loadings of the factor only
  expect_equal(sum(PTc$op == "<~" & PTc$free == 0 & PTc$ustart == 1), 2L)   # marker weights
  # Variances and intercepts of composites are computed by lavaan: fixed without a value.
  comprow <- PTc$lhs %in% c("c1", "c2") & ((PTc$op == "~~" & PTc$lhs == PTc$rhs) | PTc$op == "~1")
  expect_equal(sum(comprow), 4L)
  expect_true(all(PTc$free[comprow] == 0) && all(is.na(PTc$ustart[comprow])))
  # No covariances between indicators of different composites:
  expect_false(any(PTc$op == "~~" & PTc$lhs %in% vars9[1:3] & PTc$rhs %in% vars9[4:6]))
  # As many free parameters as the psychonetrics model:
  expect_equal(max(PTc$free), max(cvm_mod@parameters$par))
}

## Fast: guards.
expect_error(tolavaan(cnm(HS, vars = vars9, lambda = Lc, weights = Wc, beta = Bc, latent = "cov", sigma_zeta = Szc)),
             "composite = 'cov'")
expect_error(tolavaan(cvm(HS, vars = vars9, lambda = Lc, weights = Wc, latent = "ggm")), "latent = 'cov'")

## Numerical agreement with lavaan (lavaan >= 0.7 models "<~" as composites).
if (at_home() && utils::packageVersion("lavaan") >= "0.7.0"){
  lavchi <- function(l) as.numeric(lavaan::fitMeasures(l, "chisq"))
  lavdf  <- function(l) as.numeric(lavaan::fitMeasures(l, "df"))
  # Estimates of a lavaan fit in the order of a psychonetrics parameter table:
  lavest <- function(p, l, what = "est"){
    pe <- lavaan::parameterEstimates(l)
    P <- p@parameters
    lop <- c(lambda = "=~", weights = "<~", beta = "~")[P$matrix]
    idx <- match(paste(ifelse(P$matrix == "beta", P$var1, P$var2), lop, ifelse(P$matrix == "beta", P$var2, P$var1), P$group_id),
                 paste(pe$lhs, pe$op, pe$rhs, if (is.null(pe$group)) 1 else pe$group))
    pe[[what]][idx]
  }
  sel <- function(p) p@parameters$matrix %in% c("lambda", "weights", "beta") & !p@parameters$fixed

  cmod <- ' c1 <~ x1 + x2 + x3
            c2 <~ x4 + x5 + x6
            f =~ x7 + x8 + x9
            f ~ c1 + c2 '

  # (a) composite indicator covariances estimated in lavaan: everything matches.
  la <- lavaan::sem(cmod, HS, meanstructure = TRUE, likelihood = "normal", composites.cov = "free")
  pa <- fromlavaan(la, run = TRUE)
  expect_equal(pa@model, "cvm")
  expect_equal(pa@fitmeasures$df, lavdf(la))
  expect_equal(pa@fitmeasures$npar, as.numeric(lavaan::fitMeasures(la, "npar")))
  expect_true(abs(pa@fitmeasures$chisq - lavchi(la)) < 1e-4)
  expect_true(max(abs(pa@parameters$est[sel(pa)] - lavest(pa, la)[sel(pa)])) < 1e-3)
  expect_true(max(abs(pa@parameters$se[sel(pa)] - lavest(pa, la, "se")[sel(pa)])) < 1e-3)
  # The variances and intercepts of composites are computed, never estimated:
  Pa <- pa@parameters
  comp_computed <- (Pa$matrix == "sigma_zeta" & Pa$row == Pa$col & Pa$var1 %in% c("c1", "c2")) |
    (Pa$matrix == "nu_eta" & Pa$var1 %in% c("c1", "c2"))
  expect_true(all(Pa$fixed[comp_computed]))

  # (b) the lavaan default fixes these covariances to their sample values:
  #     a warning, and the same chi-square, df, estimates and standard errors
  #     of weights, loadings and regressions.
  lb <- lavaan::sem(cmod, HS, meanstructure = TRUE, likelihood = "normal")
  expect_warning(pb <- fromlavaan(lb, run = TRUE), "composites.cov")
  expect_equal(pb@fitmeasures$df, lavdf(lb))
  expect_true(abs(pb@fitmeasures$chisq - lavchi(lb)) < 1e-4)
  expect_true(max(abs(pb@parameters$est[sel(pb)] - lavest(pb, lb)[sel(pb)])) < 1e-3)
  expect_true(max(abs(pb@parameters$se[sel(pb)] - lavest(pb, lb, "se")[sel(pb)])) < 1e-3)

  # (c) a chain of composites; the intercept of the last composite equals
  #     (I - B) times the composite means (lavaan <= 0.7-3 reports another value).
  lc <- lavaan::sem(' c1 <~ x1 + x2 + x3
                      c2 <~ x4 + x5 + x6
                      c3 <~ x7 + x8 + x9
                      c2 ~ c1
                      c3 ~ c2 ', HS, meanstructure = TRUE, likelihood = "normal", composites.cov = "free")
  pc <- fromlavaan(lc, run = TRUE)
  expect_equal(pc@fitmeasures$df, lavdf(lc))
  expect_true(abs(pc@fitmeasures$chisq - lavchi(lc)) < 1e-4)
  Wpc <- as.matrix(getmatrix(pc, "weights")); mupc <- as.numeric(getmatrix(pc, "mu"))
  expect_true(max(abs(as.numeric(getmatrix(pc, "mu_eta")) - as.numeric(t(Wpc) %*% mupc))) < 1e-8)
  expect_true(max(abs(as.numeric(getmatrix(pc, "nu_eta")) -
                        as.numeric((diag(3) - as.matrix(getmatrix(pc, "beta"))) %*% t(Wpc) %*% mupc))) < 1e-8)

  # (d) two groups with equality constraints on weights (labels) and loadings.
  ld <- lavaan::sem(' c1 <~ x1 + c(a, a)*x2 + c(b, b)*x3
                      c2 <~ x4 + x5 + x6
                      f =~ x7 + x8 + x9
                      f ~ c1 + c2 ', HS, group = "school", meanstructure = TRUE, likelihood = "normal",
                    composites.cov = "free", group.equal = "loadings")
  pd <- fromlavaan(ld, run = TRUE)
  expect_equal(pd@fitmeasures$df, lavdf(ld))
  expect_true(abs(pd@fitmeasures$chisq - lavchi(ld)) < 1e-4)

  # (e) missing data.
  HSm <- HS[, vars9]
  set.seed(1); HSm[cbind(sample(nrow(HSm), 60), sample(9, 60, TRUE))] <- NA
  le <- lavaan::sem(cmod, HSm, meanstructure = TRUE, missing = "ml", composites.cov = "free", information = "expected")
  pe_ <- fromlavaan(le, run = TRUE)
  expect_equal(pe_@estimator, "FIML")
  expect_equal(pe_@fitmeasures$df, lavdf(le))
  expect_true(abs(pe_@fitmeasures$chisq - lavchi(le)) < 1e-3)

  # (f) another marker than the first indicator, from an unfitted lavaan
  #     object (psychonetrics starting values in the lavaan metric).
  mf <- ' c1 <~ NA*x1 + x2 + 1*x3
          c2 <~ NA*x4 + x5 + 2*x6
          f =~ x7 + x8 + x9
          f ~ c1 + c2 '
  pf <- fromlavaan(lavaan::sem(mf, HS, meanstructure = TRUE, likelihood = "normal", composites.cov = "free", do.fit = FALSE),
                   run = TRUE)
  Wf <- pf@parameters[pf@parameters$matrix == "weights" & pf@parameters$fixed & pf@parameters$est != 0, ]
  expect_equal(Wf$var1, c("x3", "x6"))
  expect_equal(Wf$est, c(1, 2))
  expect_equal(pf@optim$convergence, 0)
  # Equivalent to model (a). This fit starts from the psychonetrics starting
  # values, so the tolerance reflects the precision of the optimizer (the fit
  # function is flat in the direction of the composite weights):
  expect_true(abs(pf@fitmeasures$chisq - lavchi(la)) < 5e-2)

  # (g) refused models.
  expect_error(suppressWarnings(fromlavaan(lavaan::sem(sub("x1", "NA*x1", cmod), HS, composites.cov = "free"))),
               "no fixed weight")
  expect_error(suppressWarnings(fromlavaan(lavaan::sem(paste(cmod, "\n x1 ~~ x4"), HS, composites.cov = "free"))),
               "different composites")

  # (h) cvm -> lavaan: summary statistics, stored raw data, and the round trip.
  ph <- runmodel(cvm_mod)
  lh <- tolavaan(ph)
  expect_equal(ph@fitmeasures$df, lavdf(lh))
  # (tolerances reflect the precision of the psychonetrics optimizer, see (f))
  expect_true(abs(ph@fitmeasures$chisq - lavchi(lh)) < 5e-2)
  expect_true(max(abs(ph@parameters$est[sel(ph)] - lavest(ph, lh)[sel(ph)])) < 1e-2)
  expect_true(max(abs(ph@parameters$se[sel(ph)] - lavest(ph, lh, "se")[sel(ph)])) < 1e-2)
  lh2 <- tolavaan(runmodel(cvm(HS, vars = vars9, lambda = Lc, weights = Wc, beta = Bc, sigma_zeta = Szc,
                               latents = c("c1", "c2", "f"), storedata = TRUE)))
  expect_true(abs(ph@fitmeasures$chisq - lavchi(lh2)) < 5e-2)
  rt <- fromlavaan(lh, run = TRUE)
  expect_equal(rt@fitmeasures$df, ph@fitmeasures$df)
  expect_equal(rt@fitmeasures$npar, ph@fitmeasures$npar)
  expect_true(abs(rt@fitmeasures$chisq - ph@fitmeasures$chisq) < 5e-2)
}
