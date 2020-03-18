logLikelihood_Ising_group <- function(fimldata,...){
  if (missing(fimldata)){
    return(logLikelihood_Ising_group_sumstat(...))
  } else {
    return(logLikelihood_Ising_group_fiml(fimldata=fimldata,...))
  }
}

# FIML version:
logLikelihood_Ising_group_fiml <- function(...){
  stop("FIML not yet implemented for Ising distribution")
  # # Subgroup models:
  # if (!cpp){
  #   1/fulln * Reduce("+", lapply(fimldata,logLikelihood_Ising_subgroup_fiml,sigma=sigma,kappa=kappa,mu=mu))    
  # } else {
  #   1/fulln * logLikelihood_Ising_subgroup_fiml_cpp_outer(fimldata=fimldata,sigma=sigma,kappa=kappa,mu=mu)
  # }

}

# CPP version
logLikelihood_Ising_subgroup_fiml_cpp_outer <- function(...){
  stop("FIML not yet implemented for Ising distribution")
  # logLikelihood_Ising_subgroup_fiml_cpp(sigma=as.matrix(sigma),mu=as.matrix(mu),kappa = as.matrix(kappa),fimldata = fimldata,epsilon = .Machine$double.eps) 
}


logLikelihood_Ising_subgroup_fiml <- function(...){
  stop("FIML not yet implemented for Ising distribution")
  # obs <- !as.vector(dat$pattern)
  # 
  # sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
  # kappa_p <- solve_symmetric(sig_p, logdet = TRUE)
  # 
  # 
  # 
  # # Means:
  # mu_p <- mu[obs,]
  # dat$n *  logLikelihood_Ising_group_sumstat(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p)
}


# Fit function per group:
logLikelihood_Ising_group_sumstat <- function(omega,tau,beta,squares,means,responses,nobs,...){
  # Graph and thresholds:
  thresholds <- as.vector(tau)
  graph <- as.matrix(omega)
  beta <- as.numeric(beta)
  
  # Compute Z:
  Z <- computeZ(graph, thresholds, as.numeric(beta), responses)
  
  # Compute summary statistics:
  # FIXME: Not nice, will make things double
  v1 <- as.vector(means * nobs)
  v2 <- as.matrix(squares)

  
  H <-  (
    - sum(thresholds * v1) - 
      sum((graph * v2)[lower.tri(v2,diag=FALSE)])
  )
  
  # Fml
  Fml <- 2 * log(Z) + 2 * beta * H / nobs
  
  # LL:
  LL <- -(1/2)*nobs*Fml
  
  # Return:
  return(LL)
}

# Fit function for Gauss ML: -2n* log likelihood
logLikelihood_Ising <- function(prep){
  # # Prepare
  # prep <- prepareModel(parVector(model), model)
  # 
  # Add cpp:
  # for (i in seq_along(prep$groupModels)){
  #   prep$groupModels[[i]]$cpp <- model@cpp
  # }
  # 
  # Fit function per group:
  ll_per_Group <- sapply(prep$groupModels,do.call,what=logLikelihood_Ising_group)

  # Sum and return:
  sum(ll_per_Group)
}