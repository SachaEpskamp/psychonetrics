# Implied model for precision. Requires appropriate model matrices:
implied_lvm <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }
  
  if (model@cpp){
    x <- impliedcovstructures_cpp(x,"zeta",type = model@types$latent, all = all)
    x <- impliedcovstructures_cpp(x,"epsilon",type = model@types$residual, all = all)
  } else {
    x <- impliedcovstructures(x,"zeta",type = model@types$latent, all = all)
    x <- impliedcovstructures(x,"epsilon",type = model@types$residual, all = all)
  }


    # For each group:
  nGroup <- length(x)
  
  for (g in seq_along(x)){
   

    
    # Matrices I need in every model framework when estimating.
    # When beta == 0 (the default), BetaStar = (I - 0)^{-1} = I exactly and
    # t(BetaStar) %x% BetaStar = I, so the solve() and the nLat^2 x nLat^2
    # Kronecker product can be skipped (bit-identical shortcut):
      nLat_g <- nrow(x[[g]]$beta)
      if (all(x[[g]]$beta == 0)){
        BetaStar <- diag(nLat_g)
        tBetakronBeta <- diag(nLat_g^2)
      } else {
        BetaStar <- as.matrix(solve(Diagonal(nLat_g) - x[[g]]$beta))
        tBetakronBeta <- t(BetaStar) %x% BetaStar
      }
      Lambda_BetaStar <- x[[g]]$lambda %*%  BetaStar
      Betasta_sigmaZeta <- BetaStar %*% x[[g]]$sigma_zeta
      
      
      # If not all, store these extra matrices too:
      if (!all){
        x[[g]]$BetaStar <- BetaStar
        x[[g]]$Lambda_BetaStar <- Lambda_BetaStar
        x[[g]]$Betasta_sigmaZeta <- Betasta_sigmaZeta
        x[[g]]$tBetakronBeta <- tBetakronBeta
      } 
      
      # Implied means (only if mean structure is modeled):
      if (!is.null(x[[g]]$nu)){
        x[[g]]$mu <- x[[g]]$nu +  x[[g]]$lambda %*% BetaStar  %*% x[[g]]$nu_eta
      }
      
      # Implied variances:
      # Factor part:
      factorPart <- Lambda_BetaStar %*% x[[g]]$sigma_zeta %*% t(Lambda_BetaStar)

      # When corinput=TRUE, derive diagonal of sigma_epsilon from diag(sigma)=1 constraint:
      if (model@sample@corinput){
        diag(x[[g]]$sigma_epsilon) <- 1 - diag(factorPart)
      }

      x[[g]]$sigma <- factorPart + x[[g]]$sigma_epsilon

      # FIXME: forcing symmetric, but not sure why this is needed...
      x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
      
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
      
      # FIXME: forcing symmetric, but not sure why this is needed...
      # x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
      
      # # Kappa, sigma and mu never sparse:
      # x[[g]]$mu <- as.matrix(x[[g]]$mu)
      # x[[g]]$kappa <- as.matrix(x[[g]]$kappa)
      # x[[g]]$sigma <- as.matrix(x[[g]]$sigma)
  }
  
  return(x)
}
