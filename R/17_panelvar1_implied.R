# Implied model for precision. Requires appropriate model matrices:
implied_panelvar1 <- function(model,all = FALSE){
  x <- formModelMatrices(model)
  
  # For each group:
  nGroup <- length(x)
  for (g in seq_len(nGroup)){
    
    # Form the models:
    # Contemporaneous:
    if (model@types$contemporaneous == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(solve_symmetric(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"Matrix")
      }
    } else if(model@types$contemporaneous == "chol"){
      # form cov matrix:
      x[[g]]$sigma_zeta <- as(x[[g]]$lowertri_zeta %*% t(x[[g]]$lowertri_zeta), "Matrix")
      
      # Return precision and network if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(solve_symmetric(x[[g]]$sigma_zeta), "Matrix")
        x[[g]]$omega_zeta  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"Matrix")
      }
    } else if (model@types$contemporaneous == "ggm"){
      x[[g]]$sigma_zeta <- x[[g]]$delta_zeta %*% solve_symmetric(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta) %*% x[[g]]$delta_zeta
      
      # Stuff needed if all = TRUE:
      if (all){
        x[[g]]$kappa_zeta <- as(solve_symmetric(x[[g]]$sigma_zeta) , "Matrix")
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]]$IminOinv_zeta <- solve_symmetric(Diagonal(ncol(x[[g]]$omega_zeta)) - x[[g]]$omega_zeta)
      }
    } else if (model@types$contemporaneous == "prec"){
      # Precision matrix
      x[[g]]$sigma_zeta <- as(solve_symmetric(x[[g]]$kappa_zeta),"Matrix")
      
      if (all) {
        x[[g]]$omega_zeta <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_zeta)),"Matrix")
      }
    }
    
    ## Between-cases
    if (model@types$between == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        x[[g]]$kappa_mu <- as(solve_symmetric(x[[g]]$sigma_mu), "Matrix")
        x[[g]]$omega_mu  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_mu)),"Matrix")
      }
    } else if(model@types$between == "chol"){
      # form cov matrix:
      x[[g]]$sigma_mu <- as(x[[g]]$lowertri_mu %*% t(x[[g]]$lowertri_mu), "Matrix")
      
      # Return precision and network if all = TRUE:
      if (all){
        x[[g]]$kappa_mu <- as(solve_symmetric(x[[g]]$sigma_mu), "Matrix")
        x[[g]]$omega_mu  <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_mu)),"Matrix")
      }
    } else if (model@types$between == "ggm"){
      x[[g]]$sigma_mu <- x[[g]]$delta_mu %*% solve_symmetric(spectralshift(Diagonal(ncol(x[[g]]$omega_mu)) - x[[g]]$omega_mu)) %*% x[[g]]$delta_mu
      
      # Stuff needed if all = TRUE:
      if (all){
        x[[g]]$kappa_mu <- solve_symmetric(x[[g]]$sigma_mu) 
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]]$IminOinv_mu <- solve_symmetric(spectralshift(Diagonal(ncol(x[[g]]$omega_mu)) - x[[g]]$omega_mu))
      }
    } else if (model@types$between == "prec"){
      # Precision matrix
      x[[g]]$sigma_mu <- as(solve_symmetric(spectralshift(x[[g]]$kappa_mu)),"Matrix")
      
      if (all) {
        x[[g]]$omega_mu <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa_mu)),"Matrix")
      }
    }

    # Some stuff needed now:
    design <- model@extramatrices$design
    nNode <- nrow(design)
    nTime <- ncol(design)
                  
    # Identity matrix:
    I <- model@extramatrices$In
    
    # Beta star:
    # BetaStar <- as(solve_symmetric(I %x% I - (x[[g]]$beta %x% x[[g]]$beta)),"Matrix")
    BetaStar <- as(solve(I %x% I - (x[[g]]$beta %x% x[[g]]$beta)),"Matrix")
    
        
    # Vector w:
    w <- Vec(x[[g]]$sigma_mu - x[[g]]$beta %*% x[[g]]$sigma_mu %*% t(x[[g]]$beta) + x[[g]]$sigma_zeta)

    # (I - B)Sigma_mu:
    IminB_sigmu <-  (I - x[[g]]$beta) %*% x[[g]]$sigma_mu
    
    # Implied mean vector:
    fullMu <- as(do.call(rbind,lapply(seq_len(nTime),function(t){
      as.matrix(x[[g]]$mu_y)[design[,t]==1,,drop=FALSE]
    })), "Matrix")
    
    # List of implied varcovs:
    allSigmas <- list()
    allSigmas[[1]] <- Matrix(as.vector(BetaStar %*% w), nNode, nNode)

    # Fill implied:
    if (nTime > 1){
      for (t in 2:nTime){
        allSigmas[[t]] <- IminB_sigmu + x[[g]]$beta %*% allSigmas[[t-1]]
      }      
    }

    # Create the block Toeplitz:
    fullSigma  <- as(blockToeplitz(lapply(allSigmas,as.matrix)), "Matrix")
    
    # Subset and add to the list:
    x[[g]]$mu <- fullMu
    x[[g]]$sigma <- fullSigma[as.vector(design)==1,as.vector(design)==1]
  
    # FIXME: forcing symmetric, but not sure why this is needed...
    x[[g]]$sigma <- 0.5*(x[[g]]$sigma + t(x[[g]]$sigma))
    
    # Precision:
    x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    
    # FIXME: forcing symmetric, but not sure why this is needed...
    # x[[g]]$kappa <- 0.5*(x[[g]]$kappa + t(x[[g]]$kappa))
    
    
    # Implied variance--covariance:
    # sigma <- as(solve_symmetric(kappa),"dpoMatrix")
    # sigma <- as(solve_symmetric(kappa),"dpoMatrix")
    # sigma <- solve_symmetric(kappa)

    # Extra matrices needed in optimization:
    if (!all){
      x[[g]]$BetaStar <- BetaStar
      x[[g]]$E <- Emat(nrow(x[[g]]$beta),x[[g]]$beta)
      x[[g]]$w <- w
      x[[g]]$allSigmas <- allSigmas
      x[[g]]$B_Sigma_mu <- x[[g]]$beta %*% x[[g]]$sigma_mu
      x[[g]]$IkronBeta <- model@extramatrices$In %x% x[[g]]$beta
      x[[g]]$Sigma_mu_kron_I <- x[[g]]$sigma_mu %x% model@extramatrices$In
    }
    
  }

  
  x
}
