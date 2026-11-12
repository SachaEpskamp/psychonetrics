# Derivatives for the composite variable model (cvm).
#
# The Jacobian is built by the chain rule on top of the lvm derivatives. With
# the completed loading matrix lambda = lambda_f + Q (Q = T W D^-1 the composite
# loadings), V = BetaStar sigma_zeta BetaStar' and X = lambda V - T W:
#
#   d vech(Sigma) = J_lambda d vec(lambda_f)        (lvm, factor loadings)
#                 + J_Q d vec(Q)                    (lvm form, evaluated at X)
#                 + J_beta d vec(beta)              (lvm)
#                 + J_psi d vech(sigma_zeta)        (lvm)
#                 + d vech(T) - vech(Q dD Q')       (residual part)
#
# The (residual) variances / scalings of the composites are not parameters but
# follow from the constraint diag_C(V) = d = diag(W' T W). Differentiating that
# constraint (implicit function theorem) gives, for the computed quantities u,
#
#   H_C du = dd_C - R_B d vec(beta) - H_free d(free latent parameters),
#
# which is how beta, the weights, T and the free latent parameters obtain an
# additional term (-K ..., with K = J_psi U H_C^-1).

# Chain from vech(sigma) of a typed covariance block to its parameters:
cvm_typechain <- function(type, suffix, dots, L, D, Dstar, In, C_chol, A){
  n <- nrow(In)
  if (type == "cov"){
    return(Diagonal(n * (n + 1) / 2))
  } else if (type == "chol"){
    return(d_sigma_cholesky(lowertri = dots[[paste0("lowertri_", suffix)]], L = L, C = C_chol, In = In))
  } else if (type == "prec"){
    return(d_sigma_kappa(L = L, D = D, sigma = dots[[paste0("sigma_", suffix)]]))
  } else if (type == "ggm"){
    dIO <- dots[[paste0("delta_IminOinv_", suffix)]]
    return(cbind(
      d_sigma_omega(L = L, delta_IminOinv = dIO, A = A, delta = dots[[paste0("delta_", suffix)]], Dstar = Dstar),
      d_sigma_delta(L = L, delta_IminOinv = dIO, In = In, A = A, delta = dots[[paste0("delta_", suffix)]])
    ))
  } else if (type == "cor"){
    return(cbind(
      d_sigma_rho(L = L, SD = dots[[paste0("SD_", suffix)]], A = A, Dstar = Dstar),
      d_sigma_SD(L = L, SD_IplusRho = dots[[paste0("SD_IplusRho_", suffix)]], In = In, A = A)
    ))
  }
  stop("Unknown type.")
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
# All arguments are taken from "..." and extracted by their exact name: the
# group list holds many short names (S, M, D, mu, sigma, ...) that R would
# otherwise partially match to formal arguments.
d_phi_theta_cvm_group <- function(...){
  dots <- list(...)

  # Model matrices and types:
  lambda <- dots[["lambda"]]
  lambda_f <- dots[["lambda_f"]]
  weights <- dots[["weights"]]
  latent <- dots[["latent"]]
  residual <- dots[["residual"]]
  composite <- dots[["composite"]]
  meanstructure <- isTRUE(dots[["meanstructure"]])
  Tfull <- dots[["Tfull"]]
  TW <- dots[["TW"]]
  dinv <- dots[["dinv"]]
  BetaStar <- dots[["BetaStar"]]
  Lambda_BetaStar <- dots[["Lambda_BetaStar"]]
  Betasta_sigmaZeta <- dots[["Betasta_sigmaZeta"]]
  tBetakronBeta <- dots[["tBetakronBeta"]]
  sigma_eta <- dots[["sigma_eta"]]
  compvec <- dots[["compvec"]]

  # Extra matrices:
  L <- dots[["L"]]
  In <- dots[["In"]]
  C <- dots[["C"]]
  Cbeta <- dots[["Cbeta"]]
  Inlatent <- dots[["Inlatent"]]
  Deta <- dots[["Deta"]]
  L_eta <- dots[["L_eta"]]
  Dstar_eta <- dots[["Dstar_eta"]]
  Aeta <- dots[["Aeta"]]
  E_res <- dots[["E_res"]]
  Sel_res <- dots[["Sel_res"]]
  L_res <- dots[["L_res"]]
  D_res <- dots[["D_res"]]
  Dstar_res <- dots[["Dstar_res"]]
  In_res <- dots[["In_res"]]
  C_chol_res <- dots[["C_chol_res"]]
  A_res <- dots[["A_res"]]
  E_comp <- dots[["E_comp"]]
  LEE_comp <- dots[["LEE_comp"]]
  L_comp <- dots[["L_comp"]]
  D_comp <- dots[["D_comp"]]
  Dstar_comp <- dots[["Dstar_comp"]]
  In_comp <- dots[["In_comp"]]
  C_chol_comp <- dots[["C_chol_comp"]]
  A_comp <- dots[["A_comp"]]

  # Dimensions:
  nvar <- nrow(lambda)
  nlat <- ncol(lambda)
  nres <- if (is.null(E_res)) 0 else ncol(E_res)
  ncmp <- if (is.null(E_comp)) 0 else ncol(E_comp)
  nvech <- nvar * (nvar + 1) / 2

  # Composites:
  comp <- as.vector(as.matrix(compvec)) == 1
  compIdx <- which(comp)

  lambda <- as.matrix(lambda)
  lambda_f <- as.matrix(lambda_f)
  weights <- as.matrix(weights)
  BetaStar <- as.matrix(BetaStar)
  Ilat <- diag(nlat)

  # Number of rows:
  nMean <- meanstructure * nvar
  nobs <- nMean + nvech

  # Number of columns:
  nelement <- meanstructure * nvar + # nu
    meanstructure * nlat + # nu_eta
    nvar * nlat + # lambda
    nvar * nlat + # weights
    nlat^2 + # beta
    nlat * (nlat + 1) / 2 + # latent block
    nres * (nres + 1) / 2 + # residual block
    ncmp * (ncmp + 1) / 2 # composite indicator block

  # Empty Jacobian:
  Jac <- matrix(0, nobs, nelement)

  # Row indices:
  meanInds <- seq_len(nMean)
  sigmaInds <- nMean + seq_len(nvech)

  # Column indices:
  curInd <- 0
  if (meanstructure){
    nuInds <- curInd + seq_len(nvar)
    curInd <- curInd + nvar
    nuetaInds <- curInd + seq_len(nlat)
    curInd <- curInd + nlat
  }
  lambdaInds <- curInd + seq_len(nvar * nlat)
  curInd <- curInd + nvar * nlat
  weightsInds <- curInd + seq_len(nvar * nlat)
  curInd <- curInd + nvar * nlat
  betaInds <- curInd + seq_len(nlat^2)
  curInd <- curInd + nlat^2
  latentInds <- curInd + seq_len(nlat * (nlat + 1) / 2)
  curInd <- curInd + nlat * (nlat + 1) / 2
  resInds <- curInd + seq_len(nres * (nres + 1) / 2)
  curInd <- curInd + nres * (nres + 1) / 2
  compInds <- curInd + seq_len(ncmp * (ncmp + 1) / 2)

  ### Mean part ###
  if (meanstructure){
    nu <- as.matrix(dots[["nu"]])
    mu_eta <- as.matrix(dots[["mu_eta"]])
    Af <- as.matrix(dots[["Af"]])
    Mc <- as.matrix(dots[["Mc"]])

    Jac[meanInds, nuInds] <- diag(nvar) + lambda_f %*% Mc %*% t(weights)
    Jac[meanInds, nuetaInds] <- lambda_f %*% Af
    Jac[meanInds, lambdaInds] <- t(mu_eta) %x% diag(nvar)
    Jac[meanInds, weightsInds] <- lambda_f %*% Mc %*% (Ilat %x% t(nu))
    Jac[meanInds, betaInds] <- lambda_f %*% (t(mu_eta) %x% Af)
  }

  ### Variance part ###
  # lvm derivatives, evaluated at the completed matrices:
  J_lambda <- as.matrix(d_sigma_lambda_lvm(L = L, Lambda_BetaStar = Lambda_BetaStar, Betasta_sigmaZeta = Betasta_sigmaZeta, In = In, C = C))
  J_beta <- as.matrix(d_sigma_beta_lvm(L = L, lambda = lambda, Betasta_sigmaZeta = Betasta_sigmaZeta, Cbeta = Cbeta, Inlatent = Inlatent, tBetakronBeta = tBetakronBeta))
  J_psi <- as.matrix(d_sigma_sigma_zeta_lvm(L = L, Lambda_BetaStar = Lambda_BetaStar, Deta = Deta))

  # Chain to the latent parameterization:
  if (latent == "cov"){
    Jtype <- diag(nlat * (nlat + 1) / 2)
  } else if (latent == "ggm"){
    Jtype <- as.matrix(cbind(
      d_sigma_omega(L = L_eta, delta_IminOinv = dots[["delta_IminOinv_zeta"]], A = Aeta, delta = dots[["delta_zeta"]], Dstar = Dstar_eta),
      d_sigma_delta(L = L_eta, delta_IminOinv = dots[["delta_IminOinv_zeta"]], In = Inlatent, A = Aeta, delta = dots[["delta_zeta"]])
    ))
  } else stop("Latent type not supported in the cvm framework.")
  J_latent <- J_psi %*% Jtype

  if (length(compIdx) > 0){
    V <- as.matrix(sigma_eta)
    TW <- as.matrix(TW)
    Tfull <- as.matrix(Tfull)
    dv <- as.vector(as.matrix(dinv))
    Q <- lambda - lambda_f
    X <- lambda %*% V - TW

    # Position of the diagonal elements in vech(sigma_zeta):
    vechInd <- matrix(0, nlat, nlat)
    vechInd[lower.tri(vechInd, diag = TRUE)] <- seq_len(nlat * (nlat + 1) / 2)
    diagPos <- diag(vechInd)

    # U: derivative of vech(sigma_zeta) with respect to the computed quantities
    # (the variance for latent = "cov" and for composites that are isolated in
    # the latent network, the scaling delta otherwise):
    U <- matrix(0, nlat * (nlat + 1) / 2, length(compIdx))
    if (latent == "cov"){
      U[cbind(diagPos[compIdx], seq_along(compIdx))] <- 1
    } else {
      omega <- as.matrix(dots[["omega_zeta"]])
      isolated <- rowSums(omega != 0) == 0
      Jdelta <- Jtype[, nlat * (nlat - 1) / 2 + seq_len(nlat), drop = FALSE]
      for (p in seq_along(compIdx)){
        if (isolated[compIdx[p]]){
          U[diagPos[compIdx[p]], p] <- 1
        } else {
          U[, p] <- Jdelta[, compIdx[p]]
        }
      }
    }

    # Rv: derivative of diag_C(V) with respect to vech(sigma_zeta):
    Rv <- matrix(0, length(compIdx), nlat * (nlat + 1) / 2)
    # RB: derivative of diag_C(V) with respect to vec(beta):
    RB <- matrix(0, length(compIdx), nlat^2)
    for (p in seq_along(compIdx)){
      a <- BetaStar[compIdx[p], ]
      Rv[p, ] <- as.vector(t(Deta) %*% as.vector(tcrossprod(a)))
      RB[p, ] <- 2 * as.vector(outer(a, V[, compIdx[p]]))
    }

    # Implicit derivative of the computed quantities:
    HC <- Rv %*% U
    K <- (J_psi %*% U) %*% tryCatch(solve(HC), error = function(e) MASS::ginv(HC))

    J_latent <- J_latent - K %*% (Rv %*% Jtype)
    J_beta <- J_beta - K %*% RB

    # Derivative with respect to the composite loadings Q (lvm form at X):
    J_Q <- as.matrix(d_sigma_lambda_lvm(L = L, Lambda_BetaStar = X, Betasta_sigmaZeta = Ilat, In = In, C = C))

    # Derivative with respect to the composite variances d:
    J_d <- matrix(0, nvech, nlat)
    for (p in seq_along(compIdx)){
      k <- compIdx[p]
      q <- Q[, k]
      xk <- X[, k]
      M1 <- tcrossprod(q) + dv[k] * (tcrossprod(q, xk) + tcrossprod(xk, q))
      J_d[, k] <- K[, p] - M1[lower.tri(M1, diag = TRUE)]
    }

    # Weights:
    KdW <- 2 * as.matrix(t(Aeta) %*% (Ilat %x% t(TW)))
    J_W <- J_Q %*% (diag(dv, nlat, nlat) %x% Tfull) + J_d %*% KdW

    # Composite indicator covariances (T):
    Ec <- as.matrix(E_comp)
    Wc <- t(Ec) %*% weights
    J_T <- as.matrix(
      (LEE_comp + J_Q %*% ((dv * t(Wc)) %x% Ec) + J_d %*% as.matrix(t(Aeta) %*% (t(Wc) %x% t(Wc)))) %*% D_comp
    )

    # Chain to the composite parameterization:
    J_T <- J_T %*% as.matrix(cvm_typechain(composite, "c", dots, L = L_comp, D = D_comp, Dstar = Dstar_comp,
                                            In = In_comp, C_chol = C_chol_comp, A = A_comp))

    Jac[sigmaInds, weightsInds] <- J_W
    Jac[sigmaInds, compInds] <- J_T
  }

  Jac[sigmaInds, lambdaInds] <- J_lambda
  Jac[sigmaInds, betaInds] <- J_beta
  Jac[sigmaInds, latentInds] <- J_latent

  # Residuals of the factor indicators:
  if (nres > 0){
    Jac[sigmaInds, resInds] <- as.matrix(Sel_res %*% cvm_typechain(residual, "epsilon", dots, L = L_res, D = D_res, Dstar = Dstar_res,
                                                                    In = In_res, C_chol = C_chol_res, A = A_res))
  }

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_cvm <- function(prep){
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels, do.call, what = d_phi_theta_cvm_group)

  if (length(d_per_group) == 1L){
    d_per_group[[1]]
  } else {
    as.matrix(bdiag(d_per_group))
  }
}
