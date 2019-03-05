jacobian_fiml_gaussian_group_sigmaVersion_meanPart <- function(sigma,mu,kappa,data,...){
  n <- nrow(data)
  
  # Spectral shift kappa (will do nothing if kappa is pos def):
  kappa <- spectralshift(kappa)

  # Current jacobian:
  curJac <- Matrix(0,nrow=1,ncol=length(mu))
  # For every subject:
  for (i in seq_len(n)){
    if (!any(is.na(data[i,]))){
      kappa_p <- kappa
      mu_p <- mu
      y <- unlist(data[i,])
      
      # Elimination matrix:
      L <- Diagonal(n=nrow(sigma))
    } else {
      obs <- unlist(!is.na(data[i,]))
      
      # Skip if nothing to do:
      if (!any(obs)) {
        next
      }
      
      sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
      kappa_p <- corpcor::pseudoinverse(sig_p)
      
      # Handle possible non positive definiteness:
      kappa_p <- spectralshift(kappa_p)
      
      # Means:
      mu_p <- mu[obs,]
      
      # raw scores:
      y <- unlist(data[i,obs])
      
      # Find the proper elimination matrix:
      inds <- which(obs)
      L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(sigma)))
    }

    curJac <- curJac + t(y - mu_p) %*% kappa_p %*% L
    
  }
  # Mean part:
   -(2/n) * curJac
}

jacobian_fiml_gaussian_group_sigmaVersion_sigmaPart <- function(mu,sigma,D,kappa,data,...){
  n <- nrow(data)
  
  # DUmmy sigma for indices:
  dumSig <- matrix(1:(ncol(sigma)^2),ncol(sigma),ncol(sigma))
  
  # Current jacobian:
  curJac <- Matrix(0,nrow=1,ncol=nrow(sigma) * (nrow(sigma)+1) / 2)
  
  # Kappa kronecker kappa for full cases:
  KkronK <- (kappa %(x)% kappa)
  
  # For every subject:
  for (i in seq_len(n)){
    
    if (!any(is.na(data[i,]))){
      mu_p <- mu
      y <- unlist(data[i,])
      sig_p <- sigma
      # Elimination matrix:
      L <- Diagonal(n=nrow(sigma)^2)
      KkronK_p <- KkronK
    } else {
      obs <- unlist(!is.na(data[i,]))
      
      # Skip if nothing to do:
      if (!any(obs)) {
        next
      }
      
      sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
      kappa_p <- corpcor::pseudoinverse(sig_p)
      
      # Handle possible non positive definiteness:
      kappa_p <- spectralshift(kappa_p)
      
      # Means:
      mu_p <- mu[obs,]
      
      # raw scores:
      y <- unlist(data[i,obs])
      
      KkronK_p <-  (kappa_p %(x)% kappa_p)
      
      # Find the proper elimination matrix:
      inds <- c(dumSig[obs,obs,drop=FALSE])
      L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(sigma)^2))

    }
    # Update Jacobian:
    curJac <- curJac + t(Vec((y - mu_p) %*% t(y - mu_p) - sig_p)) %*% (KkronK_p) %*% L %*% D
  }
  
  # Return:
  -(1/n) * curJac
}


# jacobian function per group:
jacobian_fiml_gaussian_group_sigma <- function(...,mu,sigma){

  # Mean part:
  grad_mean <- jacobian_fiml_gaussian_group_sigmaVersion_meanPart(mu=mu,sigma=sigma,...)

  # sigma part:
  grad_sigma <- jacobian_fiml_gaussian_group_sigmaVersion_sigmaPart(mu=mu,sigma=sigma,...)

  # Combine and return:
  cbind(grad_mean,grad_sigma)
}

# Now for all groups:
jacobian_fiml_gaussian_sigma <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
