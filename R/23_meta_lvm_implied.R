# Implied model for meta-analytic LVM:
implied_meta_lvm <- function(model, all = FALSE){

  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)
  }

  # Implied covariance structures for LVM parts:
  if (model@cpp){
    x <- impliedcovstructures_cpp(x,"zeta",type = model@types$latent, all = all)
    x <- impliedcovstructures_cpp(x,"epsilon",type = model@types$residual, all = all)
    x <- impliedcovstructures_cpp(x,type = model@types$randomEffects, name = "randomEffects", all = all)
  } else {
    x <- impliedcovstructures(x,"zeta",type = model@types$latent, all = all)
    x <- impliedcovstructures(x,"epsilon",type = model@types$residual, all = all)
    x <- impliedcovstructures(x,type = model@types$randomEffects, name = "randomEffects", all = all)
  }

  for (g in seq_along(x)){

    est <- model@extramatrices$Vestimation

    # LVM implied covariance (from lvm_implied):
    BetaStar <- as.matrix(solve(Diagonal(nrow(x[[g]]$beta)) - x[[g]]$beta))
    Lambda_BetaStar <- x[[g]]$lambda %*% BetaStar
    Betasta_sigmaZeta <- BetaStar %*% x[[g]]$sigma_zeta
    tBetakronBeta <- t(BetaStar) %x% BetaStar

    # Factor part:
    factorPart <- Lambda_BetaStar %*% x[[g]]$sigma_zeta %*% t(Lambda_BetaStar)

    # Implied sigma_y (sigma_epsilon diagonal is free, so sigma_y is not necessarily a correlation matrix):
    sigma_y <- factorPart + x[[g]]$sigma_epsilon

    # Force symmetric:
    sigma_y <- 0.5 * (sigma_y + t(sigma_y))

    # Store extra matrices for derivatives:
    if (!all){
      x[[g]]$BetaStar <- BetaStar
      x[[g]]$Lambda_BetaStar <- Lambda_BetaStar
      x[[g]]$Betasta_sigmaZeta <- Betasta_sigmaZeta
      x[[g]]$tBetakronBeta <- tBetakronBeta
      x[[g]]$sigma_y <- sigma_y
    }

    # Meta-analytic structure:
    if (est == "averaged"){
      # mu = vechs(sigma_y):
      x[[g]]$mu <- Vech(sigma_y, FALSE)

      # sigma = sigma_randomEffects + V:
      x[[g]]$sigma <- x[[g]]$sigma_randomEffects + model@extramatrices[['V']]
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    } else {
      nStudy <- model@sample@groups$nobs[g]

      # Per study:
      x[[g]]$mu <- lapply(seq_len(nStudy),function(s)as.vector(Vech(sigma_y, FALSE)))

      # sigma per study:
      x[[g]]$sigma <- lapply(seq_len(nStudy),function(i) x[[g]]$sigma_randomEffects + model@extramatrices$Vall[[i]])
      x[[g]]$sigma <- lapply(x[[g]]$sigma, as.matrix)
      x[[g]]$kappa <- lapply(x[[g]]$sigma,solve_symmetric,logdet=TRUE)
      x[[g]]$kappa <- lapply(x[[g]]$kappa, as.matrix)
    }
  }

  x
}
