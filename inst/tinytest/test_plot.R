# Regression tests for plot.psychonetrics(): the qgraph wrapper added in
# 0.17.4. Everything is run with DoNotPlot = TRUE so that no device is opened;
# the returned qgraph objects are inspected instead.

suppressMessages(library(psychonetrics))

set.seed(1)
n <- 250

## ---- Data -----------------------------------------------------------------
omegaTrue <- matrix(0, 5, 5)
for (i in 1:4) omegaTrue[i, i + 1] <- omegaTrue[i + 1, i] <- 0.35
Sigma <- cov2cor(solve(diag(5) - omegaTrue))
dcs <- as.data.frame(MASS::mvrnorm(n, rep(0, 5), Sigma))
names(dcs) <- paste0("V", 1:5)

nV <- 3; Ttot <- 400
Bt <- diag(0.3, nV); Bt[2, 1] <- 0.2
ts1 <- matrix(0, Ttot, nV)
for (t in 2:Ttot) ts1[t, ] <- Bt %*% ts1[t - 1, ] + rnorm(nV)
dts <- as.data.frame(ts1); names(dts) <- paste0("X", 1:nV)

## ---- Fits -----------------------------------------------------------------
fit_ggm  <- suppressWarnings(runmodel(ggm(dcs), verbose = FALSE))
fit_cov  <- suppressWarnings(runmodel(varcov(dcs, type = "cov"), verbose = FALSE))
fit_prec <- suppressWarnings(runmodel(precision(dcs), verbose = FALSE))
fit_chol <- suppressWarnings(runmodel(cholesky(dcs), verbose = FALSE))
fit_gvar <- suppressWarnings(runmodel(gvar(dts), verbose = FALSE))

P <- function(...) suppressMessages(plot(..., DoNotPlot = TRUE))

## ---- (a) varcov defaults ---------------------------------------------------
g <- P(fit_ggm)
expect_inherits(g, "qgraph")
expect_equal(length(g$graphAttributes$Nodes$labels), 5L)
expect_equal(unname(g$graphAttributes$Nodes$labels), paste0("V", 1:5))
# bootnet-style defaults:
expect_equal(g$graphAttributes$Graph$cut, 0)
# omega is undirected and has no self-loops:
expect_false(any(g$Edgelist$directed))
expect_false(any(g$Edgelist$from == g$Edgelist$to))

## ---- (b) sigma: correlations by default, covariances on request -----------
gc1 <- P(fit_cov)
expect_true(all(gc1$Edgelist$directed))
expect_true(all(gc1$Edgelist$bidirectional))
expect_true(max(abs(gc1$Edgelist$weight)) <= 1 + 1e-8)
expect_false(any(gc1$Edgelist$from == gc1$Edgelist$to))
gc2 <- P(fit_cov, covariances = TRUE)
expect_true(any(gc2$Edgelist$from == gc2$Edgelist$to))

## ---- (c) kappa: partial correlations by default ---------------------------
gp <- P(fit_prec)
expect_false(any(gp$Edgelist$directed))
expect_true(max(abs(gp$Edgelist$weight)) <= 1 + 1e-8)
expect_true(any(P(fit_prec, covariances = TRUE)$Edgelist$from ==
                  P(fit_prec, covariances = TRUE)$Edgelist$to))

## ---- (d) Matrices without a network interpretation are refused ------------
expect_error(P(fit_chol))                 # Cholesky
expect_error(P(fit_ggm, "delta"))         # diagonal scaling
expect_error(P(fit_ggm, "mu"))            # mean vector
expect_error(P(fit_ggm, "notamatrix"))    # unknown
expect_error(P(fit_ggm, "omega", time = 1))  # 'time' needs an RI-CLPM

## ---- (e) var1 / gvar: two panels, shared layout ---------------------------
gv <- P(fit_gvar)
expect_true(is.list(gv))
expect_equal(length(gv), 2L)
expect_equal(names(gv), c("beta", "omega_zeta"))
# Temporal: directed, with the autoregressions as self-loops:
expect_true(all(gv[[1]]$Edgelist$directed))
expect_true(any(gv[[1]]$Edgelist$from == gv[[1]]$Edgelist$to))
# Contemporaneous: undirected, no self-loops:
expect_false(any(gv[[2]]$Edgelist$directed))
expect_false(any(gv[[2]]$Edgelist$from == gv[[2]]$Edgelist$to))
# Both panels are laid out identically:
expect_equal(gv[[1]]$layout, gv[[2]]$layout)
expect_equal(unname(vapply(gv, function(q) q$Arguments$title, character(1))),
             c("Temporal", "Contemporaneous"))
# Labels are the model's variable names in both panels:
expect_equal(unname(gv[[1]]$graphAttributes$Nodes$labels), paste0("X", 1:nV))
expect_equal(unname(gv[[2]]$graphAttributes$Nodes$labels), paste0("X", 1:nV))

## ---- (f) 'matrix' takes a vector and respects its order -------------------
go <- P(fit_gvar, c("omega_zeta", "beta"))
expect_equal(names(go), c("omega_zeta", "beta"))

## ---- (g) qgraph arguments override the defaults ---------------------------
gz <- P(fit_ggm, layout = "circle", cut = 0.2, labels = letters[1:5])
expect_equal(unname(gz$graphAttributes$Nodes$labels), letters[1:5])
expect_equal(gz$graphAttributes$Graph$cut, 0.2)

## ---- (h) Factor loadings --------------------------------------------------
lam <- matrix(0, 6, 2); lam[1:3, 1] <- 1; lam[4:6, 2] <- 1
Fs  <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2, 2))
dcfa <- as.data.frame(Fs %*% t(lam) + matrix(rnorm(n * 6, sd = 0.7), n, 6))
names(dcfa) <- paste0("y", 1:6)
fit_cfa <- suppressWarnings(runmodel(lvm(dcfa, lambda = lam,
                                         latents = c("F1", "F2")), verbose = FALSE))
gl <- P(fit_cfa, "lambda")
expect_inherits(gl, "qgraph")
expect_equal(unname(gl$graphAttributes$Nodes$labels), c(paste0("y", 1:6), "F1", "F2"))

## ---- (i) RI-CLPM blocks ---------------------------------------------------
pV <- 3; pT <- 3; nP <- 300
RI <- matrix(rnorm(nP * pV), nP, pV)
w  <- array(0, c(nP, pV, pT))
w[, , 1] <- matrix(rnorm(nP * pV), nP, pV)
for (t in 2:pT) w[, , t] <- w[, , t - 1] %*% t(Bt) + matrix(rnorm(nP * pV, sd = .8), nP, pV)
dpan <- do.call(cbind, lapply(1:pT, function(t) w[, , t] + RI))
colnames(dpan) <- as.vector(outer(paste0("P", 1:pV), 1:pT, paste, sep = "_"))
dpan <- as.data.frame(dpan)
des <- matrix(colnames(dpan), pV, pT); rownames(des) <- paste0("P", 1:pV)
fit_ri <- suppressWarnings(runmodel(ri_clpn(dpan, des), verbose = FALSE))
gri <- P(fit_ri)
expect_true(is.list(gri))
# Every panel is an nVar x nVar block, not the full latent-sized matrix:
expect_true(all(vapply(gri, function(q) length(q$graphAttributes$Nodes$labels),
                       numeric(1)) == pV))
expect_equal(unname(gri[[1]]$graphAttributes$Nodes$labels), paste0("P", 1:pV))
# A single wave can be selected:
expect_inherits(P(fit_ri, "omega_zeta_within", time = 2), "qgraph")
expect_true(all(P(fit_ri, "PDC", time = 1)$Edgelist$directed))
