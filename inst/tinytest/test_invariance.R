# Tests derived from psychonetrics.org examples 12 and 18 (measurement
# invariance on the bundled StarWars data).
#
# Note: the psychonetrics.org example 12 passes a grouping VECTOR via
# groups = ...; raw-data input requires a column name, so the tests add the
# grouping variable to the data and pass groups = "agegroup".
#
# Reference values hard-coded from psychonetrics 0.15.30-dev
# (fable-improvements worktree build, 2026-06-12).

suppressMessages(library(psychonetrics))
data("StarWars", package = "psychonetrics")

items <- StarWars[, 1:7]
items$agegroup <- ifelse(StarWars$Q12 < 30, "young", "older")
lambda <- matrix(0, 7, 3)
lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1; lambda[7, 3] <- 1
colnames(lambda) <- c("Original", "Prequel", "Sequel")

## ---- Example 12: configural / weak / strong invariance ----
mod_configural <- suppressWarnings(
  runmodel(lvm(items, vars = paste0("Q", 1:7), lambda = lambda,
               groups = "agegroup", identification = "variance"))
)
expect_inherits(mod_configural, "psychonetrics")
expect_true(abs(mod_configural@fitmeasures$logl - (-2914.75537)) < 0.01)
expect_equal(mod_configural@fitmeasures$df, 24)

mod_weak <- suppressWarnings(runmodel(groupequal(mod_configural, "lambda")))
expect_equal(mod_weak@fitmeasures$df, 28)

mod_strong <- suppressWarnings(runmodel(groupequal(mod_weak, "nu")))
expect_equal(mod_strong@fitmeasures$df, 32)

cmp <- compare(configural = mod_configural, weak = mod_weak, strong = mod_strong)
expect_equal(nrow(cmp), 3)
expect_true(all(is.finite(cmp$AIC)))
expect_true(all(is.finite(cmp$BIC)))
expect_true(all(is.finite(cmp$Chisq)))
# Increasingly constrained models, increasing df:
expect_equal(cmp$DF, c(24, 28, 32))
expect_true(all(is.finite(cmp$p_value[-1])))

## ---- Example 18 (reduced): equality-free modification indices ----
invisible(capture.output(mi_free <- MIs(mod_strong, type = "free", verbose = FALSE)))
expect_true(is.data.frame(mi_free))
expect_true(nrow(mi_free) > 0)
expect_true("mi_free" %in% colnames(mi_free))
expect_true(any(is.finite(mi_free$mi_free)))

## ---- at_home: example 18 fuller pipeline (10 items, partial invariance) ----
if (at_home()) {
  Data <- StarWars
  Data$agegroup <- ifelse(Data$Q12 < 30, "young", "less young")
  obsvars <- paste0("Q", 1:10)
  latents <- c("Prequels", "Original", "Sequels")
  Lambda <- matrix(0, 10, 3)
  Lambda[1:4, 1] <- 1; Lambda[c(1, 5:7), 2] <- 1; Lambda[c(1, 8:10), 3] <- 1
  Theta <- diag(1, 10); Theta[4, 10] <- Theta[10, 4] <- 1

  mod_c <- suppressWarnings(
    runmodel(lvm(Data, lambda = Lambda, vars = obsvars, latents = latents,
                 sigma_epsilon = Theta, identification = "variance",
                 groups = "agegroup"))
  )
  mod_w <- suppressWarnings(runmodel(groupequal(mod_c, "lambda")))
  mod_s <- suppressWarnings(runmodel(groupequal(mod_w, "nu")))
  expect_true(mod_c@fitmeasures$df < mod_w@fitmeasures$df)
  expect_true(mod_w@fitmeasures$df < mod_s@fitmeasures$df)

  invisible(capture.output(
    mi_nu <- MIs(mod_s, matrices = "nu", type = "free", verbose = FALSE)
  ))
  expect_true(is.data.frame(mi_nu))

  # Partial invariance: free the intercept of Q10 again:
  mod_sp <- suppressWarnings(runmodel(groupfree(mod_s, "nu", 10)))
  expect_equal(mod_sp@fitmeasures$df, mod_s@fitmeasures$df - 1)

  mod_strict <- suppressWarnings(runmodel(groupequal(mod_sp, "sigma_epsilon")))
  mod_eqvar  <- suppressWarnings(runmodel(groupequal(mod_strict, "sigma_zeta")))
  expect_true(mod_eqvar@fitmeasures$df > mod_strict@fitmeasures$df)

  cmp18 <- compare(configural = mod_c, weak = mod_w, strong = mod_s)
  expect_true(all(is.finite(cmp18$AIC)))

  nu_eta <- getmatrix(mod_eqvar, "nu_eta")
  expect_equal(length(unlist(nu_eta)), 6L)  # 3 latents x 2 groups
}
