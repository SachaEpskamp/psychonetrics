# Tests for clustersearch(): the agglomerative search for clusters of
# variables that are modeled as composites (cvm) or latent variables (lvm).
#
# The tests use population covariance matrices as input, so that the expected
# results are exact (a chi-square of zero for a true partition) rather than
# reference values from simulated data. Within a population structure several
# candidate models fit equally well, so the order in which clusters are formed
# is not tested; the tests use the clusters and the values of the criterion.

suppressMessages(library(psychonetrics))

# Shared helpers (portable multivariate normal draws; see _helpers.R):
source(if (file.exists("_helpers.R")) "_helpers.R" else
         system.file("tinytest", "_helpers.R", package = "psychonetrics"))

quiet <- function(expr) suppressMessages(expr)
asSets <- function(x) sort(sapply(x, function(v) paste(sort(v), collapse = "+")))

## ---- Population: two latent variables with three indicators each ----
lam <- matrix(0, 6, 2)
lam[1:3, 1] <- c(.8, .7, .6)
lam[4:6, 2] <- c(.75, .65, .85)
Phi <- matrix(c(1, .4, .4, 1), 2)
SigmaL <- lam %*% Phi %*% t(lam) + diag(1 - rowSums(lam^2))
# (variables on different scales, so that the input is not a correlation matrix)
sds <- c(1, 1.2, 0.9, 1.5, 0.8, 1.1)
SigmaL <- diag(sds) %*% SigmaL %*% diag(sds)
vars <- paste0("y", 1:6)
dimnames(SigmaL) <- list(vars, vars)
truth <- c("y1+y2+y3", "y4+y5+y6")

## ---- 1. Search for latent variables ----
res_l <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", cluster = "latent",
                             verbose = FALSE))
expect_inherits(res_l, "clustersearch")
expect_equal(res_l$cluster, "latent")
expect_equal(res_l$criterion, "bic")

# The search runs to a single cluster: one step per merge
expect_equal(res_l$path$step, 0:5)
expect_equal(res_l$path$nclusters, 6:1)
expect_equal(length(res_l$models), 6)
expect_equal(names(res_l$models), as.character(6:1))
expect_true(all(sapply(res_l$models, inherits, "psychonetrics")))
expect_true(all(sapply(res_l$models, function(m) m@model) == "lvm"))

# All candidate models: choose(K, 2) per step
expect_equal(as.vector(table(res_l$candidates$step)), choose(6:2, 2))
expect_equal(sum(res_l$candidates$selected), 5)
expect_true(all(res_l$candidates$eligible[res_l$candidates$selected]))

# The true partition is found, has the lowest criterion and fits exactly
expect_equal(res_l$best, 4)
expect_equal(unname(asSets(res_l$clusters)), truth)
expect_equal(unname(asSets(res_l$partitions[["2"]])), truth)
expect_true(abs(res_l$path$chisq[res_l$path$step == 4]) < 1e-3)
expect_equal(res_l$path$df[res_l$path$step == 4], 8)
expect_identical(res_l$model, res_l$models[["2"]])

# All steps up to the true partition fit exactly, the single latent variable does not
expect_true(all(abs(res_l$path$chisq[1:5]) < 1e-3))
expect_true(res_l$path$chisq[6] > 10)

# The starting model is the saturated model
expect_equal(res_l$path$chisq[1], 0)
expect_equal(res_l$path$df[1], 0)
expect_equal(res_l$path$npar[1], 27)
mod_sat <- runmodel(varcov(covs = SigmaL, means = rep(0, 6), nobs = 500, covtype = "ML"))
expect_equal(res_l$path$bic[1], mod_sat@fitmeasures$bic, tolerance = 1e-6)

# The path reports the criterion of the models
expect_equal(res_l$path$bic, unname(sapply(res_l$models, function(m) m@fitmeasures$bic)), tolerance = 1e-8)
expect_equal(res_l$path$npar, unname(sapply(res_l$models, function(m) m@fitmeasures$npar)))
expect_equal(res_l$path$chisq, unname(sapply(res_l$models, function(m) m@fitmeasures$chisq)), tolerance = 1e-3)

# Names: the cluster of step s is C<s>, other constructs carry the variable name
expect_equal(res_l$path$cluster[-1], paste0("C", 1:5))
expect_true(all(lengths(res_l$partitions[["6"]]) == 1))
expect_equal(names(res_l$partitions[["6"]]), vars)
expect_true(all(names(res_l$clusters) %in% paste0("C", 1:5)))

# Latent variables are identified by their variance by default
expect_equal(res_l$model@identification, "variance")

## ---- 2. Search for composites ----
# A latent variable model is a special case of the composite model, so the
# same structure is found. The search stops at two clusters.
res_c <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", verbose = FALSE))
expect_equal(res_c$cluster, "composite")
expect_equal(res_c$path$nclusters, 6:2)
expect_true(all(sapply(res_c$models, function(m) m@model) == "cvm"))
expect_equal(unname(asSets(res_c$partitions[["2"]])), truth)
expect_true(all(abs(res_c$path$chisq) < 1e-2))
expect_equal(as.vector(table(res_c$candidates$step)), choose(6:3, 2))
expect_true(all(res_c$candidates$eligible))
expect_true(all(!res_c$candidates$improper))
expect_equal(res_c$path$bic, unname(sapply(res_c$models, function(m) m@fitmeasures$bic)), tolerance = 1e-8)

# For two variables a composite and a latent variable are equivalent: the
# first step has the same degrees of freedom and criterion in both searches
expect_equal(res_c$path$df[2], res_l$path$df[2])
expect_equal(res_c$path$bic[2], res_l$path$bic[2], tolerance = 1e-5)
# For three variables the latent variable is more restrictive
expect_true(res_l$path$df[res_l$path$step == 4] > res_c$path$df[res_c$path$step == 4])

## ---- 3. Print and plot ----
out <- utils::capture.output(print(res_l))
expect_true(any(grepl("Cluster search", out)))
expect_true(any(grepl("latent variables", out)))
expect_true(any(grepl("Selected model: step 4", out)))
out <- utils::capture.output(print(res_c))
expect_true(any(grepl("composites", out)))

grDevices::pdf(NULL)
nodes <- plot(res_l)
expect_equal(nrow(nodes), 11)
expect_equal(sort(nodes$y[nodes$step == 0]), 1:6)
# The variables of a cluster are adjacent in the dendrogram
for (cl in res_l$clusters){
  y <- nodes$y[match(cl, nodes$name)]
  expect_equal(diff(range(y)), length(cl) - 1)
}
expect_silent(plot(res_c, labels = LETTERS[1:6], nodelabels = FALSE))
expect_error(plot(res_c, labels = "a"), "one element per variable")
grDevices::dev.off()

## ---- 4. Settings ----
# Criterion
res_aic <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", criterion = "aic",
                               min_clusters = 4, verbose = FALSE))
expect_equal(res_aic$criterion, "aic.ll")
expect_true("aic.ll" %in% names(res_aic$path))
expect_equal(res_aic$path$nclusters, 6:4)
expect_equal(res_aic$path$aic.ll, unname(sapply(res_aic$models, function(m) m@fitmeasures$aic.ll)),
             tolerance = 1e-8)

# A subset of the variables of a covariance matrix
res_v4 <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", vars = vars[c(1, 2, 4, 5)],
                              verbose = FALSE))
expect_equal(res_v4$vars, vars[c(1, 2, 4, 5)])
expect_equal(res_v4$path$nclusters, 4:2)
expect_equal(unname(asSets(res_v4$partitions[["2"]])), c("y1+y2", "y4+y5"))
expect_error(clustersearch(covs = SigmaL, nobs = 500, vars = c("y1", "y2", "z")), "names of variables in 'covs'")

# Errors
expect_error(clustersearch(covs = SigmaL, nobs = 500, lambda = diag(6)), "cannot be used in clustersearch")
expect_error(clustersearch(covs = SigmaL, nobs = 500, criterion = "rmsea"), "'criterion' must be one of")
expect_error(clustersearch(covs = SigmaL, nobs = 500, min_clusters = 1), "at least 2 in a search for composites")
expect_error(clustersearch(covs = SigmaL, nobs = 500, min_clusters = 6), "more variables")
expect_error(clustersearch(covs = SigmaL, nobs = 500, estimator = "DWLS"), "only supports the estimators")
expect_error(clustersearch(covs = SigmaL, nobs = 500, improper = NA), "'improper' must be TRUE or FALSE")
expect_error(clustersearch(nobs = 500), "Either 'data' or 'covs'")

# The function is flagged experimental with the shared note, which is shown
# once per session and is silent from version 0.19 on:
if (utils::packageVersion("psychonetrics") < "0.19"){
  exp_env <- psychonetrics:::.experimental_warned
  rm(list = grep("clustersearch", ls(exp_env), value = TRUE), envir = exp_env)
  expect_message(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", min_clusters = 5, verbose = FALSE),
                 "'clustersearch\\(\\)' is experimental in psychonetrics")
}

if (at_home()){
  ## ---- 5. Network models ----
  res_cg <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", type = "ggm", verbose = FALSE))
  expect_equal(res_cg$model@types$latent, "ggm")
  expect_equal(res_cg$model@types$composite, "ggm")
  expect_equal(res_cg$path$bic, res_c$path$bic, tolerance = 1e-6)
  expect_equal(unname(asSets(res_cg$partitions[["2"]])), truth)
  expect_true(all(sapply(res_cg$models, function(m) m@computed)))

  res_lg <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", cluster = "latent",
                                type = "ggm", verbose = FALSE))
  expect_equal(res_lg$model@types$latent, "ggm")
  expect_equal(res_lg$path$bic, res_l$path$bic, tolerance = 1e-6)
  expect_equal(unname(asSets(res_lg$clusters)), truth)

  ## ---- 6. Identification by loadings gives the same search ----
  res_ll <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", cluster = "latent",
                                identification = "loadings", verbose = FALSE))
  expect_equal(res_ll$model@identification, "loadings")
  expect_equal(res_ll$path$bic, res_l$path$bic, tolerance = 1e-5)
  expect_equal(unname(asSets(res_ll$clusters)), truth)

  ## ---- 7. Improper solutions ----
  # Candidates with an improper solution are not eligible, unless improper = TRUE
  expect_true(all(!res_l$candidates$eligible[res_l$candidates$improper]))
  res_li <- quiet(clustersearch(covs = SigmaL, nobs = 500, covtype = "ML", cluster = "latent",
                                improper = TRUE, verbose = FALSE))
  expect_true(all(res_li$candidates$eligible[res_li$candidates$converged & is.na(res_li$candidates$error)]))

  ## ---- 8. Raw data: complete data, missing data, multiple cores and groups ----
  set.seed(123)
  dat <- as.data.frame(rmvnorm(400, rep(0, 6), SigmaL))
  names(dat) <- vars
  res_d <- quiet(clustersearch(dat, verbose = FALSE))
  S <- cov(dat) * (nrow(dat) - 1) / nrow(dat)
  res_s <- quiet(clustersearch(covs = S, nobs = nrow(dat), means = colMeans(dat), covtype = "ML",
                               verbose = FALSE))
  expect_equal(res_d$path, res_s$path, tolerance = 1e-5)
  expect_equal(unname(asSets(res_d$partitions[["2"]])), truth)

  # Multiple cores give the same result
  res_d2 <- quiet(clustersearch(dat, nCores = 2, verbose = FALSE))
  expect_equal(res_d2$path, res_d$path, tolerance = 1e-8)
  expect_equal(res_d2$candidates, res_d$candidates, tolerance = 1e-8)

  # Missing data: full-information maximum likelihood
  datNA <- dat
  datNA[1:20, 1] <- NA
  datNA[21:40, 5] <- NA
  res_na <- quiet(clustersearch(datNA, cluster = "latent", verbose = FALSE))
  expect_equal(res_na$model@estimator, "FIML")
  expect_equal(unname(asSets(res_na$partitions[["2"]])), truth)

  # A subset of the variables
  res_v <- quiet(clustersearch(dat, vars = vars[1:4], verbose = FALSE))
  expect_equal(res_v$vars, vars[1:4])
  expect_equal(res_v$path$nclusters, 4:2)

  # Multiple groups
  dat$group <- rep(c("a", "b"), each = 200)
  res_g <- quiet(clustersearch(dat, groupvar = "group", cluster = "latent", verbose = FALSE))
  expect_equal(res_g$vars, vars)
  expect_equal(nrow(res_g$model@sample@groups), 2)
  expect_equal(unname(asSets(res_g$partitions[["2"]])), truth)

  # Sum scores
  res_sum <- quiet(clustersearch(dat[, vars], sumscores = TRUE, verbose = FALSE))
  w <- res_sum$model@parameters
  expect_true(all(w$fixed[w$matrix == "weights"]))
}
