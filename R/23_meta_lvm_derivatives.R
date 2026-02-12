# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_meta_lvm_group <- function(lambda, latent, residual, randomEffects, metacor, cpp, ...){

  dots <- list(...)

  # Number of observed variables:
  nvar <- nrow(lambda)

  # Number of latents:
  nlat <- ncol(lambda)

  # Number of correlations (off-diagonal):
  ncor <- nvar * (nvar - 1) / 2

  # Number of modeled elements in the mean (always correlations):
  nmod <- ncor

  # Number of observations: mean part + variance part
  nobs <- nmod + # correlations (mean part of meta-analytic model)
    nmod * (nmod+1) / 2 # Random effects (variance part)

  # Mean part and variance part indices:
  meanPart <- seq_len(nmod)
  varPart <- max(meanPart) + seq_len(nmod*(nmod+1)/2)

  #### LVM parameter count ####
  # No mean structure (always corinput), no thresholds
  # Parameters: lambda (nvar*nlat) + beta (nlat^2) + sigma_zeta (nlat*(nlat+1)/2) + sigma_epsilon (nvar*(nvar+1)/2)
  nLVM <- nvar * nlat + nlat^2 + nlat*(nlat+1)/2 + nvar*(nvar+1)/2

  #### Random effects parameter count ####
  # Depends on parameterization:
  if (randomEffects == "cov"){
    nRan <- nmod * (nmod+1) / 2
  } else if (randomEffects == "chol"){
    nRan <- nmod * (nmod+1) / 2
  } else if (randomEffects == "prec"){
    nRan <- nmod * (nmod+1) / 2
  } else if (randomEffects == "ggm"){
    nRan <- nmod * (nmod-1) / 2 + nmod  # omega (off-diag) + delta (diag)
  } else if (randomEffects == "cor"){
    nRan <- nmod * (nmod-1) / 2 + nmod  # rho (off-diag) + SD (diag)
  }

  # Total parameters:
  nTotal <- nLVM + nRan

  # Empty Jacobian:
  Jac <- matrix(0, nobs, nTotal)

  # LVM parameter indices:
  lvmInds <- seq_len(nLVM)
  ranInds <- nLVM + seq_len(nRan)

  #### Fill mean part: LVM Jacobian with corinput=TRUE ####
  # Call the existing d_phi_theta_lvm_group with corinput=TRUE and meanstructure=FALSE
  # This gives us the derivatives of vechs(sigma_y) w.r.t. LVM parameters
  # Need to rename sigma_y -> sigma for d_phi_theta_lvm_group, and remove conflicting names
  lvmDots <- dots
  lvmDots$sigma <- lvmDots$sigma_y
  # Remove names that conflict with explicit args of d_phi_theta_lvm_group:
  lvmDots$corinput <- NULL
  lvmDots$meanstructure <- NULL
  lvmDots$sigma_y <- NULL
  lvmDots$randomEffects <- NULL
  lvmDots$metacor <- NULL
  lvmDots$cpp <- NULL
  lvmDots$mu <- NULL  # meta-analytic mu would conflict with LVM mu
  lvmDots$kappa <- NULL  # meta-analytic kappa would conflict
  lvmDots$tau <- NULL  # no thresholds in meta_lvm
  lvmJac <- do.call(d_phi_theta_lvm_group, c(list(lambda = lambda, latent = latent, residual = residual,
                                                    corinput = TRUE, meanstructure = FALSE), lvmDots))

  Jac[meanPart, lvmInds] <- lvmJac

  #### Fill variance part: random effects Jacobian ####
  # Same structure as meta_varcov variance part
  if (randomEffects == "cov"){
    Jac[varPart, ranInds] <- as.matrix(Diagonal(nRan))

  } else if (randomEffects == "chol"){
    if (cpp){
      Jac[varPart, ranInds] <- d_sigma_cholesky_cpp(lowertri = dots$lowertri_randomEffects, L = dots$L_c, C = dots$C_c,
                                                     In = dots$In_c)
    } else {
      Jac[varPart, ranInds] <- d_sigma_cholesky(lowertri = dots$lowertri_randomEffects, L = dots$L_c, C = dots$C_c,
                                                 In = dots$In_c)
    }

  } else if (randomEffects == "ggm"){
    netPart <- nLVM + seq_len(nmod*(nmod-1)/2)
    scalingPart <- max(netPart) + seq_len(nmod)

    if (cpp){
      Jac[varPart, netPart] <- d_sigma_omega_cpp(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                  A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c)
      Jac[varPart, scalingPart] <- d_sigma_delta_cpp(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                      In = dots$In_c, A = dots$A_c)
    } else {
      Jac[varPart, netPart] <- d_sigma_omega(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                              A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c)
      Jac[varPart, scalingPart] <- d_sigma_delta(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                  In = dots$In_c, A = dots$A_c)
    }

  } else if (randomEffects == "prec"){
    if (cpp){
      Jac[varPart, ranInds] <- d_sigma_kappa_cpp(L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects)
    } else {
      Jac[varPart, ranInds] <- d_sigma_kappa(L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects)
    }

  } else if (randomEffects == "cor"){
    corPart <- nLVM + seq_len(nmod*(nmod-1)/2)
    sdPart <- max(corPart) + seq_len(nmod)

    if (cpp){
      Jac[varPart, corPart] <- d_sigma_rho_cpp(L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c)
      Jac[varPart, sdPart] <- d_sigma_SD_cpp(L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects, In = dots$In_c,
                                              A = dots$A_c)
    } else {
      Jac[varPart, corPart] <- d_sigma_rho(L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c)
      Jac[varPart, sdPart] <- d_sigma_SD(L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects, In = dots$In_c,
                                          A = dots$A_c)
    }
  }

  # Make sparse if needed:
  Jac <- sparseordense(Jac)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_meta_lvm <- function(prep){
  # model is already prepared!

  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_meta_lvm_group)

  # Bind by column and return:
  Reduce("bdiag",d_per_group)
}
