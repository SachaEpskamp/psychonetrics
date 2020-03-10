

jacobian_fiml_gaussian_subgroup_sigma <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- sigma[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  
  
  # Means:
  mu_p <- mu[obs]
  dat$n *  jacobian_gaussian_group_sigma(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p,D=dat$D) %*% dat$L
}



# jacobian function per group
jacobian_fiml_gaussian_group_sigma <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  # Subgroup models:
  Jac <- 1/fulln * Reduce("+", lapply(fimldata,jacobian_fiml_gaussian_subgroup_sigma,sigma=sigma,kappa=kappa,mu=mu))
  
  # Cut out the rows not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Jac <- Jac[,keep,drop=FALSE]
  }
  if (!meanstructure){
    Jac <- Jac[,-(seq_len(nvar)),drop=FALSE]
  }
  
  return(Jac)
}


# Jacobian function per group, C++ version:
jacobian_fiml_gaussian_group_sigma_cpp_outer <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  # Subgroup models:
  Jac <- 1/fulln * jacobian_fiml_gaussian_subgroup_sigma_cpp(fimldata=fimldata,sigma=as.matrix(sigma),kappa=as.matrix(kappa),mu=as.matrix(mu), epsilon = .Machine$double.eps)

  # Cut out the rows not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Jac <- Jac[,keep,drop=FALSE]
  }
  if (!meanstructure){
    Jac <- Jac[, -(seq_len(nvar)),drop=FALSE]
  }
  
  return(Jac)
}

# full FIML version:
jacobian_fiml_gaussian_group_sigma_cpp_outer_fullFIML <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  
  nDat <- length(fimldata)
  
  # Some checks:
  if (!is.list(sigma)){
    sigma <- lapply(seq_len(nDat), function(x) as.matrix(sigma))
  } else {
    if (length(sigma) != nDat){
      stop("Number of 'sigma' matrices must equal number of rows in data.")
    }  
  }
  
  if (!is.list(mu)){
    mu <- lapply(seq_len(nDat), function(x) as.vector(mu))
  } else {
    if (length(mu) != nDat){
      stop("Number of 'mu' vectors must equal number of rows in data.")
    }  
  }
  
  if (!is.list(kappa)){
    kappa <- lapply(seq_len(nDat), function(x) as.matrix(kappa))
  } else {
    if (length(kappa) != nDat){
      stop("Number of 'kappa' matrices must equal number of rows in data.")
    }  
  }
  
  # Subgroup models:
  Jac <- 1/fulln * jacobian_fiml_gaussian_subgroup_sigma_cpp_fullFIML(fimldata=fimldata,sigma=sigma,kappa=kappa,mu=mu, epsilon = .Machine$double.eps)
  
  # Cut out the rows not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Jac <- Jac[,keep]
  }
  if (!meanstructure){
    Jac <- Jac[, -(seq_len(nvar))]
  }
  
  return(Jac)
}

# 
#   # Mean part:
#   grad_mean <- jacobian_fiml_gaussian_group_sigmaVersion_meanPart(mu=mu,sigma=sigma,...)
# 
#   # sigma part:
#   grad_sigma <- jacobian_fiml_gaussian_group_sigmaVersion_sigmaPart(mu=mu,sigma=sigma,...)
# 
#   # Combine and return:
#   cbind(grad_mean,grad_sigma)
# }

# Now for all groups:
jacobian_fiml_gaussian_sigma <- function(prep){
  # model is already prepared!
  # d_phi_theta per group:
  
  # Use C++?
  if (prep$cpp){
    if (prep$fullFIML){
      # Fit function per group:
      g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma_cpp_outer_fullFIML)  
      
    } else {
      # Fit function per group:
      g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma_cpp_outer)   
    }
    
  } else {
    
    if (prep$fullFIML){
      stop("Full (rowwise) FIML only supported through C++")
    } else {
      # Fit function per group:
      g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma)
    }
    
  }
  
  
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
