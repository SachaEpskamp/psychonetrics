# Derivative of mu with respect to mu:
d_mu_mu_ggm <- function(mu,...){
  Diagonal(length(mu))
}

# Derivative of vech(kappa) with respect to vechs(omega):
d_kappa_omega_ggm <- function(L,delta,Dstar,...){
 -L %*% kronecker(delta, delta) %*% Dstar
}

# Derivative of vech(kappa) with respect to diag(delta):
d_kappa_delta_ggm <- function(L,delta,omega,A,...){
  I <- Diagonal(nrow(omega))
  L %*% (kronecker(delta %*% (I - omega), I) + kronecker(I, (I - omega)%*%delta))%*%A
}

# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_ggm_group <- function(omega,...){
  # Number of variables:
  nvar <- nrow(omega)
  
  # Number of observations:
  nobs <- nvar + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # total number of elements:
  nelement <- nvar + # Means
    (nvar * (nvar-1))/2 + # edges
    nvar # scaling

  # Empty Jacobian:
  Jac <- Matrix(0, nrow = nobs, ncol=nelement)
  
  # fill mean part:
  Jac[1:nvar,1:nvar] <- d_mu_mu_ggm(...)
  
  # Fill network part:
  Jac[nvar + seq_len((nvar * (nvar+1))/2), nvar + seq_len((nvar * (nvar-1))/2)] <- d_kappa_omega_ggm(...)

  # Fill scaling part:
  indsx <- nvar  + seq_len((nvar * (nvar+1))/2)
  indsy <- nvar + ((nvar * (nvar-1))/2) + seq_len(nvar)
  Jac[indsx,indsy] <- d_kappa_delta_ggm(omega,...)
 
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_ggm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_ggm_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}

### HESSIAN ###
# H_ggm_delta_omega_of_kappa
H_ggm_delta_omega_of_kappa <- function(delta,L,Dstar,An2,In,...){
  d <- diag(delta)
  - ( L %(x)% t(Dstar) ) %*% An2 %*% ((d %(x)% In) + (In %(x)% d))
}

# H_ggm_delta_delta_of_kappa
H_ggm_delta_delta_of_kappa <- function(omega,L,A,An2,In,In2,...){
  ones <- Matrix(1,nrow = ncol(omega),ncol=1)
  (L %(x)% t(A)) %*% (
    (In2 %(x)% (In - omega) %(x)% In) %*% An2 %*% (In %(x)% ones) + 
      (In  %(x)% (In - omega)  %(x)% In2) %*% An2 %*% (ones  %(x)% In)
  )
}

# H_ggm_omega_delta_of_kappa
H_ggm_omega_delta_of_kappa <- function(omega,L,A,An2,In,In3,Dstar,delta,E,...){
  - (L %(x)% t(A)) %*% (
    (delta %(x)% In3) %*% (In %(x)% E) %*% Dstar + 
      (In3 %(x)% delta) %*% (E %(x)% In) %*% Dstar
  )
}


# Full H-star (Hessian of model part) for a group
H_ggm_star_group <- function(...){
  # Obtain all three non-zero parts:
  H_ggm_delta_omega_of_kappa <- H_ggm_delta_omega_of_kappa(...)
  H_ggm_delta_delta_of_kappa <- H_ggm_delta_delta_of_kappa(...)
  H_ggm_omega_delta_of_kappa <- H_ggm_omega_delta_of_kappa(...)
  
  # Obtain the relevant Jacobian parts:
  j_meanPart <- jacobian_gaussian_group_kappaVersion_meanPart(...) # FIXME: Not needed except for nVar
  j_kappaPart <- jacobian_gaussian_group_kappaVersion_kappaPart(...)
  
  # Some statistics:
  nVar <- ncol(j_meanPart)
  nPar <- ncol(j_meanPart) + ncol(j_kappaPart)
  nEdge <- nVar*(nVar-1)/2
  
  # Create an empty Hessian:
  Hstar <- Matrix(0, nrow=nPar, ncol=nPar)
  
  # And Is:
  Iomega <- Diagonal(nEdge)
  Idelta <- Diagonal(nVar)
  
  # Add delta_omega_of_kappa part:
  Hstar[nVar + seq_len(nEdge),nVar + nEdge + seq_len(nVar) ] <- 
    Hstar[nVar + seq_len(nEdge),nVar + nEdge + seq_len(nVar) ] + (j_kappaPart %(x)% Iomega) %*% H_ggm_delta_omega_of_kappa
  
  # Add delta_delta_of_kappa part:
  Hstar[nVar + nEdge + seq_len(nVar),nVar + nEdge + seq_len(nVar) ] <- 
    Hstar[nVar + nEdge + seq_len(nVar),nVar + nEdge + seq_len(nVar) ] + (j_kappaPart %(x)% Idelta) %*% H_ggm_delta_delta_of_kappa
  
  # Add omega_delta_of_kappa part:
  Hstar[nVar + nEdge + seq_len(nVar),nVar + seq_len(nEdge) ] <- 
    Hstar[nVar + nEdge + seq_len(nVar),nVar + seq_len(nEdge) ] + (j_kappaPart %(x)% Idelta) %*% H_ggm_omega_delta_of_kappa
  
  # Return jacobian:
  return(Hstar)
}

# Finally, make the full Hessian for a group:
H_ggm_group <- function(...){
  # I need the Gaussian hessian:
  H_ggm_gauss <- hessian_gaussian_kappa_group(...)
  # And the parameter derivatives:
  J_ggm <- d_phi_theta_ggm_group(...)
  # And Hstar:
  Hstar <- H_ggm_star_group(...)
  
  # Now add and return:
  t(J_ggm) %*% H_ggm_gauss %*% J_ggm + Hstar
}

# Now for all groups:
H_ggm <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  H_per_group <- lapply(prep$groupModels,do.call,what=H_ggm_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    H_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * H_per_group[[i]]
  }
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",H_per_group)
}
  



