
# Per group:
expected_hessian_fiml_Gaussian_subgroup <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]

  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  
  # Means:
  mu_p <- mu[obs,]
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
  Hes <- 1/fulln * expected_hessian_fiml_Gaussian_group_cpp(fimldata=fimldata,sigma=as.matrix(sigma),kappa=as.matrix(kappa),mu=as.matrix(mu), epsilon = .Machine$double.eps)
  
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
  if (prep$cpp){
    exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group_cpp_outer)
  } else {
    exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group)    
  }

  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",exph_per_group)
}
