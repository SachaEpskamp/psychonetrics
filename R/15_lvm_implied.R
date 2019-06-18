# Implied model for precision. Requires appropriate model matrices:
implied_lvm <- function(model, all = FALSE){

  x <- formModelMatrices(model)
  
  x <- impliedcovstructures(x,"zeta",type = model@types$latent, all = all)
  x <- impliedcovstructures(x,"epsilon",type = model@types$residual, all = all)

    # For each group:
  nGroup <- length(x)
  
  for (g in seq_along(x)){
   

    
    # Matrices I need in every model framework when estimating:
      BetaStar <- as(solve(Diagonal(nrow(x[[g]]$beta)) - x[[g]]$beta),"Matrix")
      Lambda_BetaStar <- x[[g]]$lambda %*%  BetaStar 
      Betasta_sigmaZeta <- BetaStar %*% x[[g]]$sigma_zeta
      tBetakronBeta <- t(BetaStar) %x% BetaStar      
      
      
      # If not all, store these extra matrices too:
      if (!all){
        x[[g]]$BetaStar <- BetaStar
        x[[g]]$Lambda_BetaStar <- Lambda_BetaStar
        x[[g]]$Betasta_sigmaZeta <- Betasta_sigmaZeta
        x[[g]]$tBetakronBeta <- tBetakronBeta
      } 
      
      # Implied means
      x[[g]]$mu <- x[[g]]$tau +  x[[g]]$lambda %*% BetaStar  %*% x[[g]]$tau_eta
      
      # Implied variances:
      x[[g]]$sigma <- Lambda_BetaStar %*% x[[g]]$sigma_zeta %*% t(Lambda_BetaStar) + x[[g]]$sigma_epsilon
      
      # FIXME: forcing symmetric, but not sure why this is needed...
      x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
      
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
      
      # FIXME: forcing symmetric, but not sure why this is needed...
      # x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
  }
  
  return(x)
}
