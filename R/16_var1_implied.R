# Implied model for precision. Requires appropriate model matrices:
implied_var1 <- function(model,all = FALSE){
  x <- formModelMatrices(model)
  
  # For each group:
  nGroup <- length(x)
  for (g in seq_len(nGroup)){
    
    # Form the models:
    if (model@types$zeta == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(corpcor::pseudoinverse(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    } else if(model@types$zeta == "chol"){
      # form cov matrix:
      x[[g]]$sigma_zeta <- as(x[[g]]$lowertri_zeta %*% t(x[[g]]$lowertri_zeta), "Matrix")
      
      # Return precision and network if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(corpcor::pseudoinverse(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    } else if (model@types$zeta == "ggm"){
      x[[g]]$sigma_zeta <- x[[g]]$delta_zeta %*% corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta)) %*% x[[g]]$delta_zeta
      
      # Stuff needed if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- corpcor::pseudoinverse(x[[g]]$sigma_zeta) 
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]]$IminOinv_zeta <- corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta))
      }
    } else if (model@types$zeta == "prec"){
      # Precision matrix
      x[[g]]$sigma_zeta <- as(corpcor::pseudoinverse(spectralshift(x[[g]]$kappa_zeta)),"sparseMatrix")
      
      if (all) {
        x[[g]]$omega_zeta <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"sparseMatrix")
      }
    }

    # Some stuff needed now:
    BetaStar <- as(corpcor::pseudoinverse(Diagonal(nrow(x[[g]]$beta)^2) - (x[[g]]$beta %(x)% x[[g]]$beta)),"sparseMatrix")
    sigmaZetaVec <- Vec(x[[g]]$sigma_zeta)
    
    # Implied exogenous covariances:
    exoCov <- x[[g]]$exo_cholesky %*% t( x[[g]]$exo_cholesky)
    
    # Implied stationary distribution (vectorized)
    vecSigma <- BetaStar %*% sigmaZetaVec
    Sigma0 <- Matrix(as.vector(vecSigma),nrow = nrow(exoCov), ncol = ncol(exoCov))
    
    # Implied lag-1:
    Sigma1 <- x[[g]]$beta %*% Sigma0
    
    # Full implied sigma:
    x[[g]]$sigma <- rbind(
      cbind(exoCov,t(Sigma1)),
      cbind(Sigma1,Sigma0)
    )
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
    
    # Precision:
    x[[g]]$kappa <- as(corpcor::pseudoinverse(x[[g]]$sigma), "Matrix")
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
    # Implied variance--covariance:
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- as(corpcor::pseudoinverse(kappa),"dpoMatrix")
    # sigma <- corpcor::pseudoinverse(kappa)
    
    # Extra matrices needed in optimization:
    if (!all){
      x[[g]]$L_betaStar <- model@extramatrices$L %*%  BetaStar 
      x[[g]]$IkronBeta <- model@extramatrices$In %(x)% x[[g]]$beta
      x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
      x[[g]]$sigmaZetaVec <- sigmaZetaVec
      x[[g]]$BetaStar <- BetaStar
    }
    
  }

  
  x
}
