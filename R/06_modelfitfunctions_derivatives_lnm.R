# Derivative of tau with respect to mu:
d_mu_tau_lnm <- function(tau,...){
  Diagonal(length(tau))
}

# derivative of latent means:
d_mu_mueta_lnm <- function(lambda,...){
  lambda
}

# Derivative of factor loadings:
d_sigma_lambda_lnm <- function(L,lambda,sigma_eta,Ineta,C,...){
  L %*% (
    ((lambda %*% sigma_eta) %(x)% Ineta) + 
      (Ineta %(x)% (lambda %*% sigma_eta))%*%C
  )
}

# Derivative of scaling matrix:
d_sigma_delta_lnm <- function(L,lambda,IminOinv,Inlatent,A,delta_eta,...){
  L %*% (lambda %(x)% lambda) %*% (
    ((delta_eta %*% IminOinv) %(x)% Inlatent) + 
      (Inlatent %(x)% (delta_eta %*% IminOinv))
  ) %*% A
}

# Derivative of network matrix:
d_sigma_omega_lnm <- function(L,lambda,IminOinv,A,delta_eta,Dstar,...){
  LD <- lambda %*% delta_eta
  L %*% (LD %(x)% LD) %*% (IminOinv %(x)% IminOinv) %*% Dstar
}

# Derivative residual variances:
d_sigma_sigmaepsilon_lnm <- function(ITheta,...){
  ITheta
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_lnm_group <- function(lambda,...){
  # Number of variables:
  nvar <- nrow(lambda)
  
  # Number of latents:
  nlat <- ncol(lambda)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nvar + # Means
    nlat + # Latent means
    nvar * nlat + # factor loadings
    nlat*(nlat + 1)/2 + # scaling and network
    nvar *( nvar + 1)/2 # Residuals

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # Indices:
  meanInds <- 1:nvar
  sigmaInds <- nvar + seq_len(nvar*(nvar+1)/2)
  
  # Indices model:
  interceptInds <- 1:nvar
  muetaInds <- nvar + seq_len(nlat)
  lambdaInds <- max(muetaInds) + seq_len(nlat*nvar)
  omegainds <- max(lambdaInds) + seq_len(nlat*(nlat-1)/2)
  deltainds <- max(omegainds) + seq_len(nlat)
  thetainds <- max(deltainds) + seq_len(nvar*(nvar+1)/2)
  
  
  # fill intercept part:
  Jac[meanInds,interceptInds] <- d_mu_tau_lnm(...)
  
  # Fill latent mean part:
  Jac[meanInds,muetaInds] <- d_mu_mueta_lnm(lambda=lambda,...)
  
  # Fill factor loading part:
  Jac[sigmaInds,lambdaInds] <- d_sigma_lambda_lnm(lambda=lambda,...)
  
  # Fill network part:
  Jac[sigmaInds,omegainds] <- d_sigma_omega_lnm(lambda=lambda,...)
  
  # Fill scaling part:
  Jac[sigmaInds,deltainds] <- d_sigma_delta_lnm(lambda=lambda,...)
  
  # Fill residual part:
  Jac[sigmaInds,thetainds] <- d_sigma_sigmaepsilon_lnm(...)

  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_lnm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_lnm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}
# 
# ### HESSIAN ###
# # H_lnm_delta_omega_of_kappa
# H_lnm_delta_omega_of_kappa <- function(delta,L,Dstar,An2,In,...){
#   d <- diag(delta)
#   - ( L %(x)% t(Dstar) ) %*% An2 %*% ((d %(x)% In) + (In %(x)% d))
# }
# 
# # H_lnm_delta_delta_of_kappa
# H_lnm_delta_delta_of_kappa <- function(omega,L,A,An2,In,In2,...){
#   ones <- Matrix(1,nrow = ncol(omega),ncol=1)
#   (L %(x)% t(A)) %*% (
#     (In2 %(x)% (In - omega) %(x)% In) %*% An2 %*% (In %(x)% ones) + 
#       (In  %(x)% (In - omega)  %(x)% In2) %*% An2 %*% (ones  %(x)% In)
#   )
# }
# 
# # H_lnm_omega_delta_of_kappa
# H_lnm_omega_delta_of_kappa <- function(omega,L,A,An2,In,In3,Dstar,delta,E,...){
#   - (L %(x)% t(A)) %*% (
#     (delta %(x)% In3) %*% (In %(x)% E) %*% Dstar + 
#       (In3 %(x)% delta) %*% (E %(x)% In) %*% Dstar
#   )
# }
# 
# 
# # Full H-star (Hessian of model part) for a group
# H_lnm_star_group <- function(...){
#   # Obtain all three non-zero parts:
#   H_lnm_delta_omega_of_kappa <- H_lnm_delta_omega_of_kappa(...)
#   H_lnm_delta_delta_of_kappa <- H_lnm_delta_delta_of_kappa(...)
#   H_lnm_omega_delta_of_kappa <- H_lnm_omega_delta_of_kappa(...)
#   
#   # Obtain the relevant Jacobian parts:
#   j_meanPart <- jacobian_gaussian_group_kappaVersion_meanPart(...) # FIXME: Not needed except for nVar
#   j_kappaPart <- jacobian_gaussian_group_kappaVersion_kappaPart(...)
#   
#   # Some statistics:
#   nVar <- ncol(j_meanPart)
#   nPar <- ncol(j_meanPart) + ncol(j_kappaPart)
#   nEdge <- nVar*(nVar-1)/2
#   
#   # Create an empty Hessian:
#   Hstar <- Matrix(0, nrow=nPar, ncol=nPar)
#   
#   # And Is:
#   Iomega <- Diagonal(nEdge)
#   Idelta <- Diagonal(nVar)
#   
#   # Add delta_omega_of_kappa part:
#   Hstar[nVar + seq_len(nEdge),nVar + nEdge + seq_len(nVar) ] <- 
#     Hstar[nVar + seq_len(nEdge),nVar + nEdge + seq_len(nVar) ] + (j_kappaPart %(x)% Iomega) %*% H_lnm_delta_omega_of_kappa
#   
#   # Add delta_delta_of_kappa part:
#   Hstar[nVar + nEdge + seq_len(nVar),nVar + nEdge + seq_len(nVar) ] <- 
#     Hstar[nVar + nEdge + seq_len(nVar),nVar + nEdge + seq_len(nVar) ] + (j_kappaPart %(x)% Idelta) %*% H_lnm_delta_delta_of_kappa
#   
#   # Add omega_delta_of_kappa part:
#   Hstar[nVar + nEdge + seq_len(nVar),nVar + seq_len(nEdge) ] <- 
#     Hstar[nVar + nEdge + seq_len(nVar),nVar + seq_len(nEdge) ] + (j_kappaPart %(x)% Idelta) %*% H_lnm_omega_delta_of_kappa
#   
#   # Return jacobian:
#   return(Hstar)
# }
# 
# # Finally, make the full Hessian for a group:
# H_lnm_group <- function(...){
#   # I need the Gaussian hessian:
#   H_lnm_gauss <- hessian_gaussian_kappa_group(...)
#   # And the parameter derivatives:
#   J_lnm <- d_phi_theta_lnm_group(...)
#   # And Hstar:
#   Hstar <- H_lnm_star_group(...)
#   
#   # Now add and return:
#   t(J_lnm) %*% H_lnm_gauss %*% J_lnm + Hstar
# }
# 
# # Now for all groups:
# H_lnm <- function(prep){
#   # model is already prepared!
#   
#   # d_phi_theta per group:
#   H_per_group <- lapply(prep$groupModels,do.call,what=H_lnm_group)
#   
#   # Weight:
#   for (i in 1:length(prep$groupModels)){
#     H_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * H_per_group[[i]]
#   }
#   # FIXME: Computationall it is nicer to make the whole object first then fill it
#   # Bind by colum and return: 
#   Reduce("bdiag",H_per_group)
# }
#   



