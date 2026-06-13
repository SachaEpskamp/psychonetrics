# Tests for the completely standardized (std.all) solution implemented by
# parameters(x, standardized = TRUE) for the lvm and varcov families.
#
# Fast deterministic checks run always; the lavaan oracle comparisons (point
# estimates and delta-method standard errors against
# lavaan::standardizedSolution(type = "std.all")) run only at_home, because the
# small remaining differences are bounded by the optimizer-tolerance gap between
# psychonetrics' and lavaan's ML fits.

suppressMessages(library(psychonetrics))

## ---- Deterministic structural checks (no lavaan needed) ----

# CFA on HolzingerSwineford with loadings identification:
hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
lambda <- matrix(0, 9, 3)
lambda[1:3, 1] <- 1; lambda[4:6, 2] <- 1; lambda[7:9, 3] <- 1
mod_cfa <- suppressWarnings(
  runmodel(lvm(hs, lambda = lambda, identification = "loadings"), verbose = FALSE)
)

# Default behaviour unchanged: parameters() without standardized = TRUE must not
# add the std/se_std columns to the printed selection nor error:
invisible(capture.output(pt_default <- parameters(mod_cfa)))
expect_true(is.data.frame(pt_default))

# standardized = TRUE returns std and se_std on the full @parameters table:
invisible(capture.output(pt_std <- parameters(mod_cfa, standardized = TRUE)))
expect_true("std" %in% names(pt_std))
expect_true("se_std" %in% names(pt_std))

# Standardized loadings lie in a sensible range and equal the manual
# lambda * sd(eta) / sd(y) computation from the implied matrices:
imp <- mod_cfa@modelmatrices[[1]]
Sigma <- as.matrix(imp$sigma); sd_y <- sqrt(diag(Sigma))
Lambda <- as.matrix(imp$lambda)
sd_eta <- sqrt(diag(as.matrix(imp$sigma_zeta)))   # no beta -> total = disturbance
lam_rows <- pt_std[pt_std$matrix == "lambda" & pt_std$par != 0, ]
for (i in seq_len(nrow(lam_rows))){
  r <- lam_rows$row[i]; cc <- lam_rows$col[i]
  manual <- Lambda[r, cc] * sd_eta[cc] / sd_y[r]
  expect_equal(lam_rows$std[i], manual, tolerance = 1e-8)
}

# Standardized residual variances = residual / total implied variance, in (0, 1):
rv <- pt_std[pt_std$matrix == "sigma_epsilon" & pt_std$row == pt_std$col, ]
manual_rv <- diag(as.matrix(imp$sigma_epsilon)) / diag(Sigma)
expect_equal(sort(rv$std), sort(manual_rv), tolerance = 1e-8)
expect_true(all(rv$std > 0 & rv$std < 1))

# Latent variances under variance identification standardize to ~1 (diagonal of
# the latent correlation): use a variance-identified fit.
mod_cfa_var <- suppressWarnings(
  runmodel(lvm(hs, lambda = lambda, identification = "variance"), verbose = FALSE)
)
invisible(capture.output(pt_var <- parameters(mod_cfa_var, standardized = TRUE)))
lat_var <- pt_var[pt_var$matrix == "sigma_zeta" & pt_var$row == pt_var$col, ]
expect_true(all(abs(lat_var$std - 1) < 1e-6))

# varcov (plain covariance) standardized -> correlation matrix (diagonal 1):
mod_vc <- suppressWarnings(runmodel(varcov(hs[, 1:4], type = "cov"), verbose = FALSE))
invisible(capture.output(pt_vc <- parameters(mod_vc, standardized = TRUE)))
sig_diag <- pt_vc[pt_vc$matrix == "sigma" & pt_vc$row == pt_vc$col, ]
expect_true(all(abs(sig_diag$std - 1) < 1e-6))
sig_off <- pt_vc[pt_vc$matrix == "sigma" & pt_vc$row != pt_vc$col, ]
expect_true(all(abs(sig_off$std) <= 1 + 1e-8))

# Unsupported family warns and leaves std NA (use a var1 model):
# (Skip if no example data; a tiny varcov correlation-input is supported, so use
# a model framework that is NOT lvm/varcov.)
# Ising over a tiny binary dataset:
set.seed(1)
bin <- as.data.frame(matrix(sample(0:1, 60, replace = TRUE), ncol = 3))
names(bin) <- paste0("b", 1:3)
mod_is <- tryCatch(suppressWarnings(runmodel(Ising(bin), verbose = FALSE)),
                   error = function(e) NULL)
if (!is.null(mod_is)){
  expect_warning(invisible(capture.output(pt_is <- parameters(mod_is, standardized = TRUE))))
  expect_true(all(is.na(pt_is$std)))
}

## ---- lavaan oracle comparisons (at_home only) ----
if (at_home()){
  if (!requireNamespace("lavaan", quietly = TRUE)) exit_file("lavaan not available")

  std_lookup <- function(pt, mat, v1, v2 = NULL, g = 1){
    sel <- pt$matrix == mat & pt$group_id == g
    if (is.null(v2)){
      sel <- sel & pt$var1 == v1
    } else {
      sel <- sel & ((pt$var1 == v1 & pt$var2 == v2) | (pt$var1 == v2 & pt$var2 == v1))
    }
    # A loading row exists for every (observed, latent) pair; the off-target
    # ones are fixed at zero. Prefer the meaningful (non-zero) entry:
    idx <- which(sel & pt$est != 0)
    if (length(idx) == 0) idx <- which(sel)
    pt$std[idx[1]]
  }

  ## (a) CFA loadings + residual variances vs lavaan std.all:
  fit_cfa <- lavaan::cfa(
    "visual =~ x1+x2+x3\ntextual =~ x4+x5+x6\nspeed =~ x7+x8+x9",
    data = hs, meanstructure = TRUE)
  ss <- lavaan::standardizedSolution(fit_cfa, type = "std.all")
  invisible(capture.output(pt <- parameters(mod_cfa, standardized = TRUE)))
  # Loadings (psychonetrics latents are Eta_1..3 in the order visual/textual/speed):
  lav_lat <- c(visual = "Eta_1", textual = "Eta_2", speed = "Eta_3")
  for (i in seq_len(nrow(ss))){
    if (ss$op[i] == "=~"){
      pnv <- std_lookup(pt, "lambda", ss$rhs[i])
      # The std loading must match to within the optimizer-tolerance gap:
      expect_true(abs(pnv - ss$est.std[i]) < 5e-4)
    }
    if (ss$op[i] == "~~" && ss$lhs[i] == ss$rhs[i] && grepl("^x", ss$lhs[i])){
      pnv <- std_lookup(pt, "sigma_epsilon", ss$lhs[i], ss$lhs[i])
      expect_true(abs(pnv - ss$est.std[i]) < 5e-4)
    }
  }

  ## (b) SEM with beta (PoliticalDemocracy): regressions + residual covariances.
  vars <- c(paste0("y", 1:8), paste0("x", 1:3))
  pd <- lavaan::PoliticalDemocracy[, vars]
  lambda2 <- matrix(0, 11, 3)
  lambda2[9:11, 1] <- 1; lambda2[1:4, 2] <- 1; lambda2[5:8, 3] <- 1
  colnames(lambda2) <- c("ind60", "dem60", "dem65"); rownames(lambda2) <- vars
  beta <- matrix(0, 3, 3); beta[2, 1] <- 1; beta[3, 1] <- 1; beta[3, 2] <- 1
  rownames(beta) <- colnames(beta) <- c("ind60", "dem60", "dem65")
  se2 <- matrix(0, 11, 11)
  rp <- rbind(c(1, 5), c(2, 4), c(2, 6), c(3, 7), c(4, 8), c(6, 8))
  for (k in seq_len(nrow(rp))){ se2[rp[k,1], rp[k,2]] <- 1; se2[rp[k,2], rp[k,1]] <- 1 }
  diag(se2) <- 1
  mod_sem <- suppressWarnings(runmodel(
    lvm(pd, lambda = lambda2, beta = beta, sigma_epsilon = se2,
        identification = "loadings", vars = vars), verbose = FALSE))
  invisible(capture.output(pts <- parameters(mod_sem, standardized = TRUE)))
  sss <- lavaan::standardizedSolution(fit <- lavaan::sem(
    "ind60 =~ x1+x2+x3
     dem60 =~ y1+y2+y3+y4
     dem65 =~ y5+y6+y7+y8
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
     y1~~y5
     y2~~y4+y6
     y3~~y7
     y4~~y8
     y6~~y8", data = pd, meanstructure = TRUE), type = "std.all")
  beta_lat <- c(ind60 = "Eta_1", dem60 = "Eta_2", dem65 = "Eta_3")
  for (i in seq_len(nrow(sss))){
    if (sss$op[i] == "~"){   # regression -> beta
      pnv <- std_lookup(pts, "beta", beta_lat[[sss$lhs[i]]], beta_lat[[sss$rhs[i]]])
      expect_true(abs(pnv - sss$est.std[i]) < 5e-4)
    }
    if (sss$op[i] == "~~" && sss$lhs[i] != sss$rhs[i] &&
        grepl("^y", sss$lhs[i]) && grepl("^y", sss$rhs[i])){
      pnv <- std_lookup(pts, "sigma_epsilon", sss$lhs[i], sss$rhs[i])
      expect_true(abs(pnv - sss$est.std[i]) < 5e-4)   # residual correlation
    }
  }

  ## (c) Standard errors via delta method (CFA) vs lavaan std.all SEs. Only the
  ## FREE loadings carry an SE in psychonetrics; the marker (fixed) loading has
  ## se_std = NA by design (a fixed parameter has no sampling variance), whereas
  ## lavaan reports a std SE for it. Skip those.
  for (i in seq_len(nrow(ss))){
    if (ss$op[i] == "=~" && is.finite(ss$se[i]) && ss$se[i] > 0){
      sel <- pt$matrix == "lambda" & pt$var1 == ss$rhs[i] & pt$est != 0 & pt$par != 0
      idx <- which(sel)
      if (length(idx) == 0) next   # marker (fixed) loading -> no se_std
      se_pn <- pt$se_std[idx[1]]
      if (is.na(se_pn)) next
      expect_true(abs(se_pn - ss$se[i]) < 1e-3)
    }
  }
}
