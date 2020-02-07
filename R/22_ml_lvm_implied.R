# Implied model for precision. Requires appropriate model matrices:
implied_ml_lvm <- function(model,all = FALSE){
  x <- formModelMatrices(model)
  
  # Implied covariance structures:
  x <- impliedcovstructures(x, "zeta_within", type = model@types$within_latent, all = all)
  x <- impliedcovstructures(x, "epsilon_within", type = model@types$within_residual, all = all)
  x <- impliedcovstructures(x, "zeta_between", type = model@types$between_latent, all = all)
  x <- impliedcovstructures(x, "epsilon_between", type = model@types$between_residual, all = all)
  
  # For each group:
  nGroup <- length(x)
  
  
  # Some stuff needed now:
  design <- model@extramatrices$design
  nVar <- nrow(design)
  nMaxInCluster <- ncol(design)
  
  # Identity matrix for latents:
  I_eta <- model@extramatrices$I_eta
  
  for (g in 1:nGroup){
    # Beta star within:
    BetaStar_within <- as(solve(Diagonal(nrow(x[[g]]$beta_within)) - x[[g]]$beta_within),"Matrix")
    BetaStar_between <- as(solve(Diagonal(nrow(x[[g]]$beta_between)) - x[[g]]$beta_between),"Matrix")
    
    Betasta_sigmaZeta_within <- BetaStar_within %*% x[[g]]$sigma_zeta_within
    Betasta_sigmaZeta_between <- BetaStar_between %*% x[[g]]$sigma_zeta_between
    
    # 
    # Implied mean vector:
    impMu <-  x[[g]]$nu +  x[[g]]$lambda %*% BetaStar_between  %*% x[[g]]$nu_eta
    
    fullMu <- as(do.call(rbind,lapply(seq_len(nMaxInCluster),function(t){
      impMu[design[,t]==1,,drop=FALSE]
    })), "Matrix")
    
    # List of implied varcovs within-subject latents:
    nLatent <- ncol(x[[g]]$lambda)
    
    # Implied within covariance:
    sigma_eta_within <- Betasta_sigmaZeta_within %*% t(BetaStar_within)
    sigma_eta_between <- Betasta_sigmaZeta_between %*% t(BetaStar_between)
    
    # Create the block Toeplitz:
    fullSigma_within_latent  <- Diagonal(nMaxInCluster) %x% sigma_eta_within
    
    # Full within-subject cov matrix:
    fullSigma_within <- (Diagonal(nMaxInCluster) %x% x[[g]]$lambda) %*% fullSigma_within_latent %*% (Diagonal(nMaxInCluster) %x% t(x[[g]]$lambda)) + (Diagonal(nMaxInCluster) %x% x[[g]]$sigma_epsilon_within)
    
    # Full between-subject cov matrix:
    fullSigma_between <- Matrix(1,nMaxInCluster,nMaxInCluster) %x%  (
      x[[g]]$lambda %*% sigma_eta_between %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_between
    )
    
    # Full implied covmat:
    fullSigma <- fullSigma_within + fullSigma_between
    
    # Subset and add to the list:
    x[[g]]$mu <- fullMu
    x[[g]]$sigma <- fullSigma[as.vector(design)==1,as.vector(design)==1]
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
    
    # if (any(is.na( x[[g]]$sigma))){
    #   browser()
    # }
    # Precision:
    x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    # x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
    
    # Let's round to make sparse if possible:
    # x[[g]]$kappa <- as(round(x[[g]]$kappa,14),"Matrix")
    
    
    # Extra matrices needed in optimization:
    if (!all){
      Lambda_BetaStar_within <- x[[g]]$lambda %*%  BetaStar_within
      Lambda_BetaStar_between <- x[[g]]$lambda %*%  BetaStar_between
      tBetakronBeta_within <- t(BetaStar_within) %x% BetaStar_within
      tBetakronBeta_between <- t(BetaStar_between) %x% BetaStar_between
      
      x[[g]]$BetaStar_within <- BetaStar_within
      x[[g]]$BetaStar_between <- BetaStar_between
      x[[g]]$Betasta_sigmaZeta_within <- Betasta_sigmaZeta_within
      x[[g]]$Betasta_sigmaZeta_between <- Betasta_sigmaZeta_between
      x[[g]]$Lambda_BetaStar_within <- Lambda_BetaStar_within
      x[[g]]$Lambda_BetaStar_between <- Lambda_BetaStar_between
      x[[g]]$tBetakronBeta_within <- tBetakronBeta_within
      x[[g]]$tBetakronBeta_between <- tBetakronBeta_between
      
      # x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
      x[[g]]$sigma_eta_within <- sigma_eta_within
      x[[g]]$sigma_eta_between <- sigma_eta_between
      x[[g]]$lamWkronlamW <- x[[g]]$lambda %x% x[[g]]$lambda
    } else {
      x[[g]]$sigma_within <- x[[g]]$lambda %*% sigma_eta_within %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_within
      x[[g]]$sigma_between <- x[[g]]$lambda %*% sigma_eta_between %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_between
      x[[g]]$sigma_eta_within <- sigma_eta_within
      x[[g]]$sigma_eta_between <- sigma_eta_between
      x[[g]]$sigma_eta <- sigma_eta_within + sigma_eta_between
      x[[g]]$sigma_crosssection <- x[[g]]$sigma_within + x[[g]]$sigma_between
      
      
      
    }
    
  }
  
  x
}
