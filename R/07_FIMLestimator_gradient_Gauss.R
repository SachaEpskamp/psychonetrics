

jacobian_fiml_gaussian_subgroup_sigma <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- sigma[obs,obs,drop=FALSE]
  kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  
  
  
  
  # Means:
  mu_p <- mu[obs,]
  dat$n *  jacobian_gaussian_group_sigma(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p,D=dat$D) %*% dat$L
}



# jacobian function per group
jacobian_fiml_gaussian_group_sigma <- function(fimldata,fulln,sigma,kappa,mu,...){
  
  # Subgroup models:
  1/fulln * Reduce("+", lapply(fimldata,jacobian_fiml_gaussian_subgroup_sigma,sigma=sigma,kappa=kappa,mu=mu))
  
}


# Jacobian function per group, C++ version:
jacobian_fiml_gaussian_group_sigma_cpp_outer <- function(fimldata,fulln,sigma,kappa,mu,...){
  # Subgroup models:
  1/fulln * jacobian_fiml_gaussian_subgroup_sigma_cpp(fimldata=fimldata,sigma=as.matrix(sigma),kappa=as.matrix(kappa),mu=as.matrix(mu), epsilon = .Machine$double.eps)
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
  if (prep$cpp){
    g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma_cpp_outer)  
  } else {
    g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_fiml_gaussian_group_sigma)
  }
  
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
