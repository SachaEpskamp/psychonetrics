# Tests for sampling weights (pseudo-ML) via the sampling_weights= argument of
# varcov()/lvm(). Sampling weights force the robust MLR configuration (Huber-White
# sandwich SEs + Yuan-Bentler-Mplus scaled test), matching lavaan's
# sampling.weights= (Asparouhov 2005; Rosseel 2012). Structural checks run always;
# the numerical agreement with lavaan runs at_home().

suppressMessages(library(psychonetrics))

## ---- structural checks (no lavaan) ----
data("HolzingerSwineford1939", package = "lavaan")
HS <- HolzingerSwineford1939
set.seed(1); HS$w <- runif(nrow(HS), 0.3, 3)
vars <- paste0("x", 1:9)
Lam <- matrix(0, 9, 3); Lam[1:3,1] <- Lam[4:6,2] <- Lam[7:9,3] <- 1

m <- runmodel(lvm(HS, lambda = Lam, vars = vars, sampling_weights = "w",
                  identification = "loadings"), verbose = FALSE)
# Sampling weights map to the internal ML estimator + MLR robust config:
expect_equal(m@estimator, "ML")
expect_equal(psychonetrics:::get_robust_config(m)$label, "MLR")
# Normalized weights are stored and (single group) sum to n:
sw <- psychonetrics:::get_sample_weights(m@sample)
expect_true(length(sw) == 1)
expect_equal(sum(sw[[1]]), nrow(HS))
expect_true(psychonetrics:::has_sampling_weights(m))

# Error guards:
expect_error(lvm(HS, lambda = Lam, vars = vars, ordered = vars, sampling_weights = "w"),
             "ordinal")
expect_error(varcov(HS, vars = vars, sampling_weights = "not_a_column"),
             "not found")
HSm <- HS; HSm$x1[1:20] <- NA
expect_error(lvm(HSm, lambda = Lam, vars = vars, sampling_weights = "w", estimator = "FIML"),
             "missing")

## ---- constant weights reduce to unweighted MLR ----
HS$w1 <- 1
mc <- runmodel(lvm(HS, lambda = Lam, vars = vars, sampling_weights = "w1",
                   identification = "loadings"), verbose = FALSE)
mu <- runmodel(lvm(HS, lambda = Lam, vars = vars, estimator = "MLR",
                   identification = "loadings"), verbose = FALSE)
pc <- mc@parameters[!mc@parameters$fixed, ]
pu <- mu@parameters[!mu@parameters$fixed, ]
expect_equal(pc$est, pu$est, tolerance = 1e-6)
expect_equal(pc$se,  pu$se,  tolerance = 1e-6)

## ---- numerical agreement with lavaan (at_home) ----
if (at_home() && requireNamespace("lavaan", quietly = TRUE)){
  suppressMessages(library(lavaan))
  library(dplyr)
  lavmod <- 'f1 =~ x1+x2+x3
f2 =~ x4+x5+x6
f3 =~ x7+x8+x9'

  freelam <- function(mod) mod@parameters %>%
    dplyr::filter(matrix == "lambda", !fixed) %>% dplyr::arrange(group_id, col, row)
  lavlam <- function(lf){
    pe <- parameterEstimates(lf)
    if (is.null(pe$group)) pe$group <- 1L
    pe %>% dplyr::filter(op == "=~", !rhs %in% c("x1","x4","x7")) %>%
      dplyr::mutate(v = as.integer(sub("x","",rhs)), f = match(lhs, c("f1","f2","f3"))) %>%
      dplyr::arrange(group, f, v)
  }

  # (1) single-group weighted CFA vs lavaan sampling.weights:
  m1 <- runmodel(lvm(HS, lambda = Lam, vars = vars, sampling_weights = "w",
                     identification = "loadings"), verbose = FALSE)
  l1 <- cfa(lavmod, HS, sampling.weights = "w", meanstructure = TRUE)
  p1 <- freelam(m1); e1 <- lavlam(l1)
  expect_equal(p1$est, e1$est, tolerance = 1e-3)
  expect_equal(p1$se,  e1$se,  tolerance = 1e-3)           # robust Huber-White SEs
  lfm <- fitMeasures(l1)
  expect_equal(m1@fitmeasures$chisq.scaled, as.numeric(lfm["chisq.scaled"]), tolerance = 1e-2)
  expect_equal(m1@fitmeasures$chisq.scaling.factor, as.numeric(lfm["chisq.scaling.factor"]), tolerance = 1e-3)
  expect_equal(m1@fitmeasures$cfi.robust, as.numeric(lfm["cfi.robust"]), tolerance = 1e-3)

  # (2) two-group weighted CFA vs lavaan (global weight normalization):
  m2 <- runmodel(lvm(HS, lambda = Lam, vars = vars, groups = "school",
                     sampling_weights = "w", identification = "loadings"), verbose = FALSE)
  l2 <- cfa(lavmod, HS, group = "school", sampling.weights = "w", meanstructure = TRUE)
  p2 <- freelam(m2); e2 <- lavlam(l2)
  expect_equal(p2$est, e2$est, tolerance = 2e-3)
  expect_equal(p2$se,  e2$se,  tolerance = 3e-3)
  expect_equal(m2@fitmeasures$chisq, as.numeric(fitMeasures(l2)["chisq"]), tolerance = 1e-2)
  # weighted per-group covariances match lavaan's sampstat:
  ss <- lavInspect(l2, "sampstat")
  expect_equal(as.numeric(m2@sample@covs[[1]]), as.numeric(ss[[1]]$cov), tolerance = 1e-8)

  # (3) weighted SEM with regressions vs lavaan:
  data("PoliticalDemocracy", package = "lavaan"); PD <- PoliticalDemocracy
  set.seed(7); PD$w <- runif(nrow(PD), 0.5, 2.5)
  Lsem <- matrix(0,11,3); Lsem[1:3,1]<-1; Lsem[4:7,2]<-1; Lsem[8:11,3]<-1
  Bsem <- matrix(0,3,3); Bsem[2,1]<-Bsem[3,1]<-Bsem[3,2]<-1
  mS <- runmodel(lvm(PD, lambda = Lsem, vars = c(paste0("x",1:3), paste0("y",1:8)),
                     beta = Bsem, latents = c("ind60","dem60","dem65"),
                     sigma_zeta = diag(3), identification = "loadings",
                     sampling_weights = "w"), verbose = FALSE)
  lS <- sem("ind60=~x1+x2+x3\ndem60=~y1+y2+y3+y4\ndem65=~y5+y6+y7+y8\ndem60~ind60\ndem65~ind60+dem60",
            PD, sampling.weights = "w", meanstructure = TRUE)
  pB <- sort(mS@parameters$est[mS@parameters$matrix=="beta" & !mS@parameters$fixed])
  eB <- sort(parameterEstimates(lS)$est[parameterEstimates(lS)$op=="~"])
  pBse <- sort(mS@parameters$se[mS@parameters$matrix=="beta" & !mS@parameters$fixed])
  eBse <- sort(parameterEstimates(lS)$se[parameterEstimates(lS)$op=="~"])
  expect_equal(pB, eB, tolerance = 2e-3)
  expect_equal(pBse, eBse, tolerance = 3e-3)
}
