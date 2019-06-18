# Implied model for precision. Requires appropriate model matrices:
implied_tsdlvm1 <- function(model,all = FALSE){
  x <- formModelMatrices(model)
  
  # Implied covariance structures:
  x <- impliedcovstructures(x, "zeta", type = model@types$zeta, all = all)
  x <- impliedcovstructures(x, "epsilon", type = model@types$epsilon, all = all)
  # x <- impliedcovstructures(x, "zeta_between", type = model@types$between_latent, all = all)
  # x <- impliedcovstructures(x, "epsilon_between", type = model@types$between_residual, all = all)
  
  # For each group:
  nGroup <- length(x)
  
  
  # Some stuff needed now:
  
  # Identity matrix for latents:
  I_eta <- model@extramatrices$I_eta
  
  for (g in 1:nGroup){
    # Beta star:
    BetaStar <- as(solve(I_eta %x% I_eta - (x[[g]]$beta %x% x[[g]]$beta)),"Matrix")
    
    # Implied mean vector:
    impMu <- x[[g]]$tau + x[[g]]$lambda %*% x[[g]]$mu_eta
    
    fullMu <- as(rbind(x[[g]]$exo_means,impMu), "Matrix")
    
    # Exogenous cov part:
    exoCov <- x[[g]]$exo_cholesky %*% t( x[[g]]$exo_cholesky)
    
    # Latent lag-0:
    nLatent <- ncol(x[[g]]$lambda)
    Sigma_eta_0 <- Matrix(as.vector(BetaStar %*% Vec(x[[g]]$sigma_zeta)), nLatent, nLatent)
    
    # Observed stationary:
    Sigma_y_0 <-  x[[g]]$lambda %*%  Sigma_eta_0 %*% t(x[[g]]$lambda) + x[[g]]$sigma_epsilon
    
    # Lag 1 part:
    Sigma_eta_1 <- x[[g]]$beta %*% Sigma_eta_0
    
    # Lag 1 observed:
    Sigma_y_1 <-  x[[g]]$lambda %*%  Sigma_eta_1 %*% t(x[[g]]$lambda) 
    
    
    # Subset and add to the list:
    x[[g]]$mu <- fullMu
    
    # Full implied sigma:
    x[[g]]$sigma <- rbind(
      cbind(exoCov,t(Sigma_y_1)),
      cbind(Sigma_y_1,Sigma_y_0)
    )
    
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
      x[[g]]$Sigma_eta_0 <- Sigma_eta_0
      x[[g]]$Sigma_eta_1 <- Sigma_eta_1
      x[[g]]$IkronBeta <- model@extramatrices$I_eta %x% x[[g]]$beta
      x[[g]]$lamWkronlamW <- x[[g]]$lambda %x% x[[g]]$lambda
    } else {
      # Add PDC:
      x[[g]]$PDC <- computePDC(x[[g]]$beta,x[[g]]$kappa_zeta)
    }
    
  }
  
  
  x
}
