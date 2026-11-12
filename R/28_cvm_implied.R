# Implied model for the composite variable model (cvm).
#
# The cvm framework implements the factor- and composite-based SEM of
# Schamberger, Schuberth, Henseler & Rosseel: observed variables are either
# indicators of common factors (y_f = Lambda eta_f + epsilon) or composite
# indicators forming a composite (eta_c = W' y_c). The model is evaluated by
# rewriting the composite part as an ordinary latent variable model:
#
#   T      = block-diagonal covariance matrix of the composite indicators (the
#            "_c" matrices: sigma_c, omega_c, ...), embedded in a P x P matrix
#   D      = diag(W' T W)                 (the composite variances)
#   Q      = T W D^-1                     (the composite loadings)
#   Theta  = T - T W D^-1 W' T            (singular residual covariance)
#
# after which Sigma = (Lambda + Q) V (Lambda + Q)' + Theta_f + Theta, with
# V = (I - B)^-1 Sigma_zeta (I - B)^-T. The diagonal of V must equal D for the
# composites, which determines the (residual) variances of the composites in
# sigma_zeta (or their scalings in delta_zeta): these are computed here and are
# not free parameters. The same holds for the composite intercepts in nu_eta.

# Solve for the computed (residual) variances / scalings of the composites.
# Returns the latent matrices with the computed entries filled in.
#
# The computed residual variance of a composite is allowed to become negative
# (the explained variance then exceeds the variance of the composite): the
# implied model is simply continued smoothly, exactly as a free variance
# parameter would be, so that the R and C++ pipelines evaluate the same
# function and its analytic derivatives remain valid everywhere. Such a
# solution is reported with a warning after estimation.
cvm_solve_composites <- function(grp, comp, d, B, BetaStar, latent){
  nLat <- length(comp)
  proper <- TRUE

  if (latent == "cov"){
    S0 <- as.matrix(grp$sigma_zeta)
    if (any(comp)){
      diag(S0)[comp] <- 0
      G <- BetaStar[comp, comp, drop = FALSE]^2
      rhs <- d[comp] - diag(BetaStar %*% S0 %*% t(BetaStar))[comp]
      psi <- tryCatch(as.vector(solve(G, rhs)), error = function(e) as.vector(MASS::ginv(G) %*% rhs))
      if (any(!is.finite(psi)) || any(psi <= 0)) proper <- FALSE
      diag(S0)[comp] <- psi
    }
    grp$sigma_zeta <- S0

  } else if (latent == "ggm"){
    omega <- as.matrix(grp$omega_zeta)
    delta <- diag(as.matrix(grp$delta_zeta))

    # Residual variances of composites that are isolated in the network; these
    # are placed directly in sigma_zeta after it is formed (NA: not used):
    varpatch <- rep(NA_real_, nLat)

    if (any(comp)){
      IminOinv <- as.matrix(solve_symmetric(diag(nLat) - omega))
      eps <- sqrt(.Machine$double.eps)

      # Composites that are isolated in the network get their residual variance
      # from the linear system; all other composites must be exogenous, in which
      # case the scaling follows in closed form from the composite variance:
      isolated <- comp & (rowSums(omega != 0) == 0)
      connected <- comp & !isolated

      delta[comp] <- 0
      if (any(connected)){
        # Exogenous: no incoming effects (tested on beta itself, which is exact):
        if (any(B[connected, , drop = FALSE] != 0)){
          stop("A composite cannot have both incoming effects (beta) and edges in the latent network (omega_zeta).")
        }
        val <- d[connected] / diag(IminOinv)[connected]
        bad <- !is.finite(val) | val <= 0
        if (any(bad)){
          # Only possible when a covariance matrix in the model is not positive definite:
          proper <- FALSE
          val[bad] <- eps
        }
        delta[connected] <- sqrt(val)
      }
      if (any(isolated)){
        S0 <- (delta * IminOinv) %*% diag(delta, nLat, nLat)
        G <- BetaStar[isolated, isolated, drop = FALSE]^2
        rhs <- d[isolated] - diag(BetaStar %*% S0 %*% t(BetaStar))[isolated]
        psi <- tryCatch(as.vector(solve(G, rhs)), error = function(e) as.vector(MASS::ginv(G) %*% rhs))
        if (any(!is.finite(psi)) || any(psi <= 0)) proper <- FALSE
        varpatch[isolated] <- psi
        # The scaling is only defined for a positive variance:
        psipos <- psi
        psipos[!is.finite(psipos) | psipos <= 0] <- eps
        delta[isolated] <- sqrt(psipos)
      }
    }
    grp$delta_zeta <- diag(delta, nLat, nLat)
    grp$compvar_patch <- matrix(varpatch, ncol = 1)
  } else {
    stop("Latent type '", latent, "' is not supported in the cvm framework.")
  }

  grp$proper_composites <- proper
  grp
}

# Core function, working on a list of formed model matrices:
implied_cvm_core <- function(x, types, extramatrices, meanstructure = TRUE, all = FALSE){
  latent <- types$latent
  residual <- types$residual
  composite <- types$composite

  # Which blocks are in the model?
  has_res <- !is.null(extramatrices$E_res)
  has_comp <- !is.null(extramatrices$E_comp)

  # Which constructs are composites?
  comp <- as.vector(as.matrix(extramatrices$compvec)) == 1

  # Implied covariance structures of the residual and composite indicator blocks:
  if (has_res){
    x <- impliedcovstructures(x, "epsilon", type = residual, all = all)
  }
  if (has_comp){
    x <- impliedcovstructures(x, "c", type = composite, all = all)
  }

  ### Step 1: composite variances and the computed latent entries ###
  for (g in seq_along(x)){
    W <- as.matrix(x[[g]]$weights)
    B <- as.matrix(x[[g]]$beta)
    nLat <- ncol(W)
    nVar <- nrow(W)

    # Composite indicator covariances (T), embedded in a nVar x nVar matrix:
    if (has_comp){
      Ec <- as.matrix(extramatrices$E_comp)
      Tfull <- Ec %*% as.matrix(x[[g]]$sigma_c) %*% t(Ec)
      TW <- Tfull %*% W
      d <- colSums(W * TW)
      d[!comp] <- 0
      dinv <- ifelse(comp, 1/d, 0)
      dinv[!is.finite(dinv)] <- 0
    } else {
      Tfull <- matrix(0, nVar, nVar)
      TW <- matrix(0, nVar, nLat)
      d <- rep(0, nLat)
      dinv <- rep(0, nLat)
    }

    # (I - B)^-1:
    if (all(B == 0)){
      BetaStar <- diag(nLat)
    } else {
      BetaStar <- as.matrix(solve(diag(nLat) - B))
    }

    # Computed variances / scalings:
    x[[g]] <- cvm_solve_composites(x[[g]], comp, d, B, BetaStar, latent)

    # Store:
    x[[g]]$Tfull <- Tfull
    x[[g]]$TW <- TW
    x[[g]]$dvec <- matrix(d, ncol = 1)
    x[[g]]$dinv <- matrix(dinv, ncol = 1)
    x[[g]]$BetaStar <- BetaStar
  }

  # Implied latent covariance structure (now that the computed entries are in):
  x <- impliedcovstructures(x, "zeta", type = latent, all = all)

  ### Step 2: implied means and covariances ###
  for (g in seq_along(x)){
    W <- as.matrix(x[[g]]$weights)
    B <- as.matrix(x[[g]]$beta)
    nLat <- ncol(W)
    nVar <- nrow(W)
    Tfull <- x[[g]]$Tfull
    TW <- x[[g]]$TW
    dinv <- as.vector(x[[g]]$dinv)
    BetaStar <- x[[g]]$BetaStar
    sigma_zeta <- as.matrix(x[[g]]$sigma_zeta)

    # latent = "ggm": the residual variance of a composite that is isolated in
    # the network is placed in sigma_zeta directly (identical to the squared
    # scaling when positive, and a smooth continuation otherwise):
    if (!is.null(x[[g]]$compvar_patch)){
      patch <- as.vector(x[[g]]$compvar_patch)
      idx <- which(!is.na(patch))
      if (length(idx) > 0){
        sigma_zeta[cbind(idx, idx)] <- patch[idx]
        x[[g]]$sigma_zeta <- sigma_zeta
      }
    }

    # Composite loadings:
    Q <- t(t(TW) * dinv)

    # Factor loadings only (composite columns are computed):
    lambda_f <- as.matrix(x[[g]]$lambda)
    lambda_f[, comp] <- 0

    # Full loading matrix:
    lambda <- lambda_f + Q

    Lambda_BetaStar <- lambda %*% BetaStar
    Betasta_sigmaZeta <- BetaStar %*% sigma_zeta
    sigma_eta <- Betasta_sigmaZeta %*% t(BetaStar)

    # Residual part:
    residPart <- Tfull - Q %*% t(TW)
    if (has_res){
      Ef <- as.matrix(extramatrices$E_res)
      residPart <- residPart + Ef %*% as.matrix(x[[g]]$sigma_epsilon) %*% t(Ef)
    }

    # Implied covariances:
    sigma <- Lambda_BetaStar %*% sigma_zeta %*% t(Lambda_BetaStar) + residPart
    sigma <- 0.5 * (sigma + t(sigma))

    # Mean structure. Composite means follow from the indicator means
    # (W' nu); the means of the factors from the structural model:
    if (!is.null(x[[g]]$nu)){
      nu <- as.matrix(x[[g]]$nu)
      nu_eta <- as.matrix(x[[g]]$nu_eta)
      Af <- matrix(0, nLat, nLat)
      if (any(!comp)){
        Af[!comp, !comp] <- solve(diag(sum(!comp)) - B[!comp, !comp, drop = FALSE])
      }
      Mc <- diag(nLat) + Af %*% B
      Mc[, !comp] <- 0
      mu_eta <- Af %*% nu_eta + Mc %*% (t(W) %*% nu)
      mu <- nu + lambda_f %*% mu_eta
    }

    # Store:
    x[[g]]$lambda <- lambda
    x[[g]]$sigma <- sigma
    x[[g]]$kappa <- solve_symmetric(sigma, logdet = TRUE)
    x[[g]]$sigma_eta <- sigma_eta
    if (!is.null(x[[g]]$nu)){
      x[[g]]$mu <- mu
      x[[g]]$mu_eta <- mu_eta
    }

    if (!all){
      # Extra matrices needed for the derivatives:
      x[[g]]$lambda_f <- lambda_f
      x[[g]]$Lambda_BetaStar <- Lambda_BetaStar
      x[[g]]$Betasta_sigmaZeta <- Betasta_sigmaZeta
      if (all(B == 0)){
        x[[g]]$tBetakronBeta <- diag(nLat^2)
      } else {
        x[[g]]$tBetakronBeta <- t(BetaStar) %x% BetaStar
      }
      if (!is.null(x[[g]]$nu)){
        x[[g]]$Af <- Af
        x[[g]]$Mc <- Mc
      }
    } else {
      # Computed composite intercepts:
      if (!is.null(x[[g]]$nu)){
        x[[g]]$nu_eta <- (diag(nLat) - B) %*% mu_eta
      }
      # Implied (singular) residual covariances of the composite indicators:
      if (has_comp){
        Ec <- as.matrix(extramatrices$E_comp)
        x[[g]]$sigma_epsilon_c <- t(Ec) %*% (Tfull - Q %*% t(TW)) %*% Ec
      }
      # Remove the internal matrices:
      x[[g]]$Tfull <- NULL
      x[[g]]$TW <- NULL
      x[[g]]$dvec <- NULL
      x[[g]]$dinv <- NULL
      x[[g]]$BetaStar <- NULL
      x[[g]]$proper_composites <- NULL
      x[[g]]$compvar_patch <- NULL
    }
  }

  return(x)
}

# Implied model for the cvm framework:
implied_cvm <- function(model, all = FALSE){
  x <- formModelMatrices(model)
  implied_cvm_core(x, types = model@types, extramatrices = model@extramatrices,
                   meanstructure = model@meanstructure, all = all)
}
