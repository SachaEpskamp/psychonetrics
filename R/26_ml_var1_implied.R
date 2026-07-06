# Implied model for ml_var1. Produces, per group, the two-level distribution
# parameters consumed by the Gauss2L estimator: mu (2p), sigma_within (2p x 2p),
# sigma_between (2p x 2p). R-only (model@cpp is forced FALSE by the constructor).
implied_ml_var1 <- function(model, all = FALSE){
  x <- formModelMatrices(model)

  toeplitz <- isTRUE(model@types$toeplitz)

  # Implied covariance structures for the two typed blocks:
  x <- impliedcovstructures(x, "zeta_within",  type = model@types$within_latent,  all = all)
  x <- impliedcovstructures(x, "zeta_between", type = model@types$between_latent, all = all)

  nGroup <- length(x)

  for (g in seq_len(nGroup)){

    if (toeplitz){
      beta <- as.matrix(x[[g]]$beta)
      p <- nrow(beta)
      In <- model@extramatrices$In

      # Stationary within structure:
      BetaStar <- as.matrix(solve(Diagonal(p^2) - (beta %x% beta)))
      sigmaZetaVec <- Vec(x[[g]]$sigma_zeta_within)
      Sigma0 <- matrix(as.vector(BetaStar %*% sigmaZetaVec), p, p)
      Sigma0 <- 0.5 * (Sigma0 + t(Sigma0))
      Sigma1 <- beta %*% Sigma0

      SigmaW <- rbind(
        cbind(Sigma0, t(Sigma1)),
        cbind(Sigma1, Sigma0)
      )
      SigmaW <- as.matrix(0.5 * (SigmaW + t(SigmaW)))

      # Two-level distribution parameters:
      mu <- as.matrix(x[[g]]$mu)
      x[[g]]$mu <- rbind(mu, mu)                                      # 2p x 1
      x[[g]]$sigma_within <- SigmaW                                   # 2p x 2p
      x[[g]]$sigma_between <- as.matrix(matrix(1, 2, 2) %x% x[[g]]$sigma_zeta_between)

      if (!all){
        x[[g]]$BetaStar <- BetaStar
        x[[g]]$IkronBeta <- In %x% beta
        x[[g]]$sigmaZetaVec <- sigmaZetaVec
      } else {
        x[[g]]$sigma0 <- Sigma0
        x[[g]]$sigma1 <- Sigma1
        # Partial directed correlations (temporal network), if kappa available:
        if (!is.null(x[[g]]$kappa_zeta_within)){
          x[[g]]$PDC <- computePDC(beta, x[[g]]$kappa_zeta_within)
        }
      }

    } else {
      # Free-within saturated mode: the 2p within/between blocks are the typed
      # matrices directly (chol at dimension 2p), mu is the free 2p mean:
      x[[g]]$mu <- as.matrix(x[[g]]$mu)
      x[[g]]$sigma_within <- as.matrix(x[[g]]$sigma_zeta_within)
      x[[g]]$sigma_between <- as.matrix(x[[g]]$sigma_zeta_between)
    }
  }

  x
}
