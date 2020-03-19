
# Per group:
expected_hessian_fiml_Gaussian_subgroup <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]

  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  
  # Means:
  mu_p <- mu[obs]
  dat$n * t(dat$L) %*% expected_hessian_Gaussian_group(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p,D=dat$D) %*% dat$L
}

# Fit function per group:
expected_hessian_fiml_Gaussian_group <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  
  # Subgroup models:
  Hes <- 1/fulln * Reduce("+", lapply(fimldata,expected_hessian_fiml_Gaussian_subgroup,sigma=sigma,kappa=kappa,mu=mu))
  
  # Cut out the rows and columns not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Hes <- Hes[keep,keep]
  }
  if (!meanstructure){
    Hes <- Hes[-(seq_len(nvar)),-(seq_len(nvar))]
  }
  
  return(Hes)
}

# C++ version
expected_hessian_fiml_Gaussian_group_cpp_outer <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  # Subgroup models:
  Hes <- 1/fulln * expected_hessian_fiml_Gaussian_group_cppversion(fimldata=fimldata,sigma=as.matrix(sigma),kappa=as.matrix(kappa),mu=as.matrix(mu), epsilon = .Machine$double.eps)
  
  # Cut out the rows and columns not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Hes <- Hes[keep,keep]
  }
  if (!meanstructure){
    Hes <- Hes[-(seq_len(nvar)),-(seq_len(nvar))]
  }
  
  return(Hes)
}

# Full FIML version:
expected_hessian_fiml_Gaussian_group_cpp_outer_fullFIML <- function(fimldata,fulln,sigma,kappa,mu,means, meanstructure = TRUE, corinput = FALSE,...){
  
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
  Hes <- 1/fulln * expected_hessian_fiml_Gaussian_group_cpp_fullFIML(fimldata=fimldata,sigma=sigma,kappa=kappa,mu=mu, epsilon = .Machine$double.eps)
  
  # Cut out the rows and columns not needed
  # FIXME: Nicer to not have to compute these in the first place...
  nvar <- ncol(sigma)
  if (corinput){
    keep <- c(rep(TRUE,nvar),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Hes <- Hes[keep,keep]
  }
  if (!meanstructure){
    Hes <- Hes[-(seq_len(nvar)),-(seq_len(nvar))]
  }
  
  return(Hes)
}

# expected_hessian_fiml_Gaussian_group <- function(...){
# 
#   # Mean part:
#   meanPart <- expected_hessian_fiml_Gaussian_group_meanPart(...)
#  
#   # Variance part:
#   varPart <- expected_hessian_fiml_Gaussian_group_varPart(...)
#   
#   # Return as block matrix:
#   foo <- bdiag(meanPart, varPart)
# }

# Total:
expected_hessian_fiml_Gaussian <- function(prep){
  # model is already prepared!
  # d_phi_theta per group:
  
  # Use C++?
  if (prep$cpp){
    if (prep$fullFIML){
      # Fit function per group:
      exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group_cpp_outer_fullFIML)   
      
    } else {
      # Fit function per group:
      exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group_cpp_outer)      
    }
    
  } else {
    
    if (prep$fullFIML){
      stop("Full (rowwise) FIML only supported through C++")
    } else {
      # Fit function per group:
      exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group)   
    }
    
  }
  
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",exph_per_group)
}
