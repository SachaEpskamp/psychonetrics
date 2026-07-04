# Implied model for panelvar. Requires appropriate model matrices.
# The panelvar model is the observed-variable special case of dlvm1
# (lambda = I, zero residual variances, observed means): the implied
# structures below are the dlvm1 ones with all factor-loading algebra
# stripped out.
implied_panelvar <- function(model,all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)
  }

  if (model@cpp){
    # Implied covariance structures:
    x <- impliedcovstructures_cpp(x, "zeta_within", type = model@types$within_latent, all = all)
    x <- impliedcovstructures_cpp(x, "zeta_between", type = model@types$between_latent, all = all)
  } else {
    # Implied covariance structures:
    x <- impliedcovstructures(x, "zeta_within", type = model@types$within_latent, all = all)
    x <- impliedcovstructures(x, "zeta_between", type = model@types$between_latent, all = all)
  }

  # For each group:
  nGroup <- length(x)

  # Some stuff needed now:
  design <- model@extramatrices$design
  nVar <- nrow(design)
  nTime <- ncol(design)

  # Identity matrix for the variables:
  In <- model@extramatrices$In

  for (g in 1:nGroup){
    # Beta star:
    BetaStar <- as.matrix(solve(In %x% In - (x[[g]]$beta %x% x[[g]]$beta)))

    # Implied mean vector: the stationary means, repeated at every wave and
    # subset by the design matrix:
    mu <- x[[g]]$mu
    fullMu <- as(do.call(rbind,lapply(seq_len(nTime),function(t){
      mu[design[,t]==1,,drop=FALSE]
    })), "Matrix")

    # List of implied within-person varcovs per lag:
    allSigmas_within <- list()
    allSigmas_within[[1]] <- matrix(as.vector(BetaStar %*% Vec(x[[g]]$sigma_zeta_within)), nVar, nVar)

    # Fill implied:
    if (nTime > 1){
      for (t in 2:nTime){
        allSigmas_within[[t]] <- x[[g]]$beta %*% allSigmas_within[[t-1]]
      }
    }

    # Create the block Toeplitz (full within-person cov matrix):
    fullSigma_within <- blockToeplitz(lapply(allSigmas_within,as.matrix))

    # Full between-person cov matrix:
    fullSigma_between <- Matrix(1,nTime,nTime) %x% x[[g]]$sigma_zeta_between

    # Full implied covmat:
    fullSigma <- fullSigma_within + fullSigma_between

    # Subset and add to the list:
    x[[g]]$mu <- as.matrix(fullMu)
    x[[g]]$sigma <- fullSigma[as.vector(design)==1,as.vector(design)==1]

    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$sigma <- as.matrix(0.5*(x[[g]]$sigma + t(x[[g]]$sigma)))

    # Precision:
    x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)

    # Extra matrices needed in optimization:
    if (!all){
      x[[g]]$BetaStar <- BetaStar
      x[[g]]$allSigmas_within <- allSigmas_within
      x[[g]]$IkronBeta <- In %x% x[[g]]$beta
    } else {
      # Named as in dlvm1 (with lambda = I the observed-level and
      # "latent"-level structures coincide) so that downstream users of
      # getmatrix() keep working:
      x[[g]]$sigma_within <- as.matrix(allSigmas_within[[1]])
      x[[g]]$sigma_between <- as.matrix(x[[g]]$sigma_zeta_between)
      x[[g]]$sigma_within_full <- as.matrix(fullSigma_within)
      x[[g]]$sigma_eta_within <- as.matrix(allSigmas_within[[1]])
      x[[g]]$sigma_eta_within_lag1 <- as.matrix(allSigmas_within[[2]])
      x[[g]]$sigma_crosssection <- as.matrix(allSigmas_within[[1]] + x[[g]]$sigma_zeta_between)

      # Add PDC:
      if (!is.null(x[[g]]$kappa_zeta_within)){
        x[[g]]$PDC <- computePDC(x[[g]]$beta,x[[g]]$kappa_zeta_within)
      } else {
        x[[g]]$PDC <- computePDC(x[[g]]$beta,solve_symmetric(x[[g]]$sigma_zeta_within))
      }

    }
  }

  x
}
