# Implied model for precision. Requires appropriate model matrices:
implied_dlvm1 <- function(model,all = FALSE){
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
  nTime <- ncol(design)
  
  # Identity matrix for latents:
  I_eta <- model@extramatrices$I_eta
  
  for (g in 1:nGroup){
    # Beta star:
    BetaStar <- as(solve(I_eta %x% I_eta - (x[[g]]$beta %x% x[[g]]$beta)),"Matrix")
    
    # Implied mean vector:
    impMu <- x[[g]]$tau + x[[g]]$lambda %*% x[[g]]$mu_eta
    
    fullMu <- as(do.call(rbind,lapply(seq_len(nTime),function(t){
      impMu[design[,t]==1,,drop=FALSE]
    })), "Matrix")
    
    # List of implied varcovs within-subject latents:
    nLatent <- ncol(x[[g]]$lambda)
    
    allSigmas_within <- list()
    allSigmas_within[[1]] <- Matrix(as.vector(BetaStar %*% Vec(x[[g]]$sigma_zeta_within)), nLatent, nLatent)
    
    # Fill implied:
    if (nTime > 1){
      for (t in 2:nTime){
        allSigmas_within[[t]] <- x[[g]]$beta %*% allSigmas_within[[t-1]]
      }      
    }
    
    # Create the block Toeplitz:
    fullSigma_within_latent  <- as(blockToeplitz(lapply(allSigmas_within,as.matrix)), "Matrix")
    
    # Full within-subject cov matrix:
    fullSigma_within <- (Diagonal(nTime) %x% x[[g]]$lambda) %*% fullSigma_within_latent %*% (Diagonal(nTime) %x% t(x[[g]]$lambda)) + (Diagonal(nTime) %x% x[[g]]$sigma_epsilon_within)
    
    # Full between-subject cov matrix:
    fullSigma_between <- Matrix(1,nTime,nTime) %x%  (
      x[[g]]$lambda %*% x[[g]]$sigma_zeta_between %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_between
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
      x[[g]]$BetaStar <- BetaStar
      # x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
      x[[g]]$allSigmas_within <- allSigmas_within
      x[[g]]$IkronBeta <- model@extramatrices$I_eta %x% x[[g]]$beta
      x[[g]]$lamWkronlamW <- x[[g]]$lambda %x% x[[g]]$lambda
    } else {
      x[[g]]$sigma_within <- x[[g]]$lambda %*% allSigmas_within[[1]] %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_between
      x[[g]]$sigma_between <- x[[g]]$lambda %*% x[[g]]$sigma_zeta_between %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon_between
      x[[g]]$sigma_within_full <- fullSigma_within

      # Add PDC:
      x[[g]]$PDC <- computePDC(x[[g]]$beta,x[[g]]$kappa_zeta_within)
    }
    
  }
  
  x
}
