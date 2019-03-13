# Implied model for precision. Requires appropriate model matrices:
implied_lvm <- function(model, all = FALSE){

  x <- formModelMatrices(model)

    # For each group:
  nGroup <- length(x)
  
  for (g in seq_along(x)){
    
    ### LATENT VARIANCES ###
    if (model@types$latent == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(corpcor::pseudoinverse(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    } else if(model@types$latent == "chol"){
      # form cov matrix:
      x[[g]]$sigma_zeta <- as(x[[g]]$lowertri_zeta %*% t(x[[g]]$lowertri_zeta), "dsyMatrix")
      
     # Return precision and network if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(corpcor::pseudoinverse(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    } else if (model@types$latent == "ggm"){
      x[[g]]$sigma_zeta <- x[[g]]$delta_zeta %*% corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta)) %*% x[[g]]$delta_zeta
      
      # Stuff needed if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- corpcor::pseudoinverse(x[[g]]$sigma_zeta) 
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]]$IminOinv_zeta <- corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta))
      }
    } else if (model@types$latent == "prec"){
      # Precision matrix
      x[[g]]$sigma_zeta <- corpcor::pseudoinverse(spectralshift(x[[g]]$kappa_zeta))
      
      if (all) {
        x[[g]]$omega_zeta <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    }
    
    ### RESIDUAL VARIANCES ###
    if (model@types$residual == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        x[[g]]$kappa_epsilon <- as(corpcor::pseudoinverse(x[[g]]$sigma_epsilon), "Matrix")
        x[[g]]$omega_epsilon  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_epsilon)),"sparseMatrix")
      }
    } else if(model@types$residual == "chol"){
      # form cov matrix:
      x[[g]]$sigma_epsilon <- as(x[[g]]$lowertri_epsilon %*% t(x[[g]]$lowertri_epsilon), "Matrix")
      
      # Return precision and network if all = TRUE:
      if (all){
        x[[g]]$kappa_epsilon <- as(corpcor::pseudoinverse(x[[g]]$sigma_epsilon), "Matrix")
        x[[g]]$omega_epsilon  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_epsilon)),"sparseMatrix")
      }
    } else if (model@types$residual == "ggm"){
      x[[g]]$sigma_epsilon <- x[[g]]$delta_epsilon %*% corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_epsilon)) - x[[g]]$omega_epsilon)) %*% x[[g]]$delta_epsilon
      
      # Stuff needed if all = TRUE:
      if (all){
        x[[g]]$kappa_epsilon <- corpcor::pseudoinverse(x[[g]]$sigma_epsilon) 
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]]$IminOinv_epsilon <- corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_epsilon)) - x[[g]]$omega_epsilon))
      }
    } else if (model@types$residual == "prec"){
      # Precision matrix
      x[[g]]$sigma_epsilon <- corpcor::pseudoinverse(spectralshift(x[[g]]$kappa_epsilon))
      
      if (all) {
        x[[g]]$omega_epsilon <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_epsilon)),"sparseMatrix")
      }
    }
    

    
    # Matrices I need in every model framework when estimating:
      BetaStar <- as(corpcor::pseudoinverse(Diagonal(nrow(x[[g]]$beta)) - x[[g]]$beta),"sparseMatrix")
      Lambda_BetaStar <- x[[g]]$lambda %*%  BetaStar 
      Betasta_sigmaZeta <- BetaStar %*% x[[g]]$sigma_zeta
      tBetakronBeta <- t(BetaStar) %(x)% BetaStar      
      
      
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
      
      x[[g]]$kappa <- corpcor::pseudoinverse(x[[g]]$sigma)
      
      # FIXME: forcing symmetric, but not sure why this is needed...
      x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
  }
  
  return(x)
}
