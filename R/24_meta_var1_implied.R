# Implied model for meta-analytic VAR(1):
implied_meta_var1 <- function(model, all = FALSE){

  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)
  }

  # Implied covariance structures for contemporaneous and random effects:
  if (model@cpp){
    x <- impliedcovstructures_cpp(x, "zeta", type = model@types$zeta, all = all)
    x <- impliedcovstructures_cpp(x, type = model@types$randomEffects, name = "randomEffects", all = all)
  } else {
    x <- impliedcovstructures(x, "zeta", type = model@types$zeta, all = all)
    x <- impliedcovstructures(x, type = model@types$randomEffects, name = "randomEffects", all = all)
  }

  for (g in seq_along(x)){

    est <- model@extramatrices$Vestimation
    nNode <- nrow(x[[g]]$beta)

    # VAR(1) implied covariance (same math as var1, but no exoCov):
    BetaStar <- as.matrix(solve(Diagonal(nNode^2) - (x[[g]]$beta %x% x[[g]]$beta)))
    sigmaZetaVec <- Vec(x[[g]]$sigma_zeta)

    # Implied stationary distribution (vectorized):
    vecSigma0 <- BetaStar %*% sigmaZetaVec
    Sigma0 <- matrix(as.vector(vecSigma0), nrow = nNode, ncol = nNode)

    # Implied lag-1 cross-covariance:
    Sigma1 <- x[[g]]$beta %*% Sigma0

    # Force symmetric on Sigma0:
    Sigma0 <- 0.5 * (Sigma0 + t(Sigma0))

    # Store intermediates for derivatives:
    if (!all){
      x[[g]]$BetaStar <- BetaStar
      x[[g]]$sigmaZetaVec <- sigmaZetaVec
      x[[g]]$IkronBeta <- Diagonal(nNode) %x% x[[g]]$beta
      x[[g]]$L_betaStar <- model@extramatrices$L %*% BetaStar
      x[[g]]$Sigma0 <- Sigma0
      x[[g]]$Sigma1 <- Sigma1
    }

    # Meta-analytic mean = [vech(Sigma0), vec(Sigma1)]:
    muVec <- c(Vech(Sigma0, diag = TRUE), as.vector(Sigma1))

    # Meta-analytic structure:
    if (est == "averaged"){
      x[[g]]$mu <- muVec

      # sigma = sigma_randomEffects + V:
      x[[g]]$sigma <- x[[g]]$sigma_randomEffects + model@extramatrices[['V']]
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    } else {
      nStudy <- model@sample@groups$nobs[g]

      # Per study:
      x[[g]]$mu <- lapply(seq_len(nStudy), function(s) muVec)

      # sigma per study:
      x[[g]]$sigma <- lapply(seq_len(nStudy), function(i) x[[g]]$sigma_randomEffects + model@extramatrices$Vall[[i]])
      x[[g]]$sigma <- lapply(x[[g]]$sigma, as.matrix)
      x[[g]]$kappa <- lapply(x[[g]]$sigma, solve_symmetric, logdet = TRUE)
      x[[g]]$kappa <- lapply(x[[g]]$kappa, as.matrix)
    }
  }

  x
}
