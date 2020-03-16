# # Cholesky derivative:
# d_sigma_cholesky <- function(lowertri,L,C,In,...){
#   
#   L %*% ((In %x% In) + C) %*% ((lowertri %x% In) %*% t(L))
# }
# 
# # Derivative of scaling matrix:
# d_sigma_delta <- function(L,delta_IminOinv,In,A,delta,...){
#   L %*% (
#     (delta_IminOinv%x% In) + 
#       (In %x% delta_IminOinv)
#   ) %*% A
# }
# 
# # Derivative of network matrix:
# d_sigma_omega <- function(L,delta_IminOinv,A,delta,Dstar,...){
#   # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
#   
#   # delta_IminOinv <- delta %*% IminOinv
#   L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
#   
#   # all(a == b)
# }
# 
# # Derivative of precision matrix:
# d_sigma_kappa <- function(L,D,sigma,...){
#   - L %*% (sigma %x% sigma) %*% D
# }
# 
# # Derivative of rho:
# d_sigma_rho <- function(L,SD,A,delta,Dstar,...){
#   L %*% (SD %x% SD) %*% Dstar
# }
# 
# # Derivative of SDs:
# d_sigma_SD <- function(L,SD_IplusRho,In,A,...){
#   L %*% (
#     (SD_IplusRho%x% In) + 
#       (In %x% SD_IplusRho)
#   ) %*% A
# }
# 
# # Derivative of omega to covariances in the metacor = TRUE setting:
# d_sigma_omega_corinput <- function(L,delta_IminOinv,A,delta,Dstar,IminOinv,In,...){
#   # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
#   
#   # delta_IminOinv <- delta %*% IminOinv
#   L %*% (
#     (delta_IminOinv %x% delta_IminOinv) -
#       1/2 *  ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A %*% 
#       Diagonal(x = diag(IminOinv)^(-1.5)) %*% t(A) %*% (IminOinv %x% IminOinv)
#   ) %*% Dstar
# }



# Full jacobian of phi (distribution parameters) with respect to theta (model parameters) for a group
d_phi_theta_meta_varcov_group <- function(y,randomEffects,metacor,cpp,...){
  
  # Dots:
  dots <- list(...)
  
  # Number of variables:
  nvar <- nrow(dots$sigma_y)
  
  # Number of correlations:
  ncor <- nvar * (nvar-1) / 2
  
  # Number of modeled elements:
  nmod <- ncor + !metacor * nvar
  
  # Number of observations:
  nobs <- nmod + # correlations
    nmod * (nmod+1) / 2 # Random effects
  
  # Mean part:
  meanPart <- seq_len(nmod)
  
  # Variance part:
  varPart <- max(meanPart) + seq_len(nmod*(nmod+1)/2)    
  
  # Empty Jacobian:
  Jac <- matrix(0, nobs, nobs)
  
  
  
  # Fill the mean part with model for cors/covs:
  if (metacor){
    Lmat <- dots$Lstar 
  } else {
    Lmat <- dots$L
  }
  
  if (y == "cov"){
    # Regular covs:
    Jac[meanPart,meanPart] <- as.matrix(Diagonal(nmod))
  } else if (y == "chol"){
    # Cholesky decomposition:
    if (cpp){
      Jac[meanPart,meanPart] <- d_sigma_cholesky_cpp(lowertri = dots$lowertri_y,L = Lmat, C = dots$C, In = dots$In )      
    } else {
      Jac[meanPart,meanPart] <- d_sigma_cholesky(lowertri = dots$lowertri_y,L = Lmat, C = dots$C, In = dots$In )      
    }
    
  } else if (y == "ggm"){
    
    # Gaussian graphical model:
    netPart <- seq_len(nvar*(nvar-1)/2)
    scalingPart <- max(netPart) + seq_len(nvar)
    
    if (metacor){
      if (cpp){
        Jac[meanPart,netPart] <- d_sigma_omega_corinput_cpp(delta_IminOinv = dots$delta_IminOinv_y, 
                                                            L = Lmat, A = dots$A, delta = dots$delta_y,
                                                            Dstar = dots$Dstar, IminOinv = dots$IminOinv_y, 
                                                            In = dots$In)        
      } else {
        Jac[meanPart,netPart] <- d_sigma_omega_corinput(delta_IminOinv = dots$delta_IminOinv_y, 
                                                        L = Lmat, A = dots$A, delta = dots$delta_y,
                                                        Dstar = dots$Dstar, IminOinv = dots$IminOinv_y, 
                                                        In = dots$In)
      }
      
      
    } else {
      if (cpp){
        Jac[meanPart,netPart] <- d_sigma_omega_cpp(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
                                                   A = dots$A, delta = dots$delta_y, Dstar = dots$Dstar)
        # Jac[meanPart,scalingPart] <- d_sigma_delta_cpp(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
        #                                                In = dots$In, A = dots$A, delta = dots$delta_y)
        Jac[meanPart,scalingPart] <- d_sigma_delta_cpp(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
                                                       In = dots$In, A = dots$A)
      } else {
        Jac[meanPart,netPart] <- d_sigma_omega(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
                                               A = dots$A, delta = dots$delta_y, Dstar = dots$Dstar)
        # Jac[meanPart,scalingPart] <- d_sigma_delta(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
        #                                            In = dots$In, A = dots$A, delta = dots$delta_y)
        Jac[meanPart,scalingPart] <- d_sigma_delta(L = Lmat, delta_IminOinv = dots$delta_IminOinv_y,
                                                   In = dots$In, A = dots$A)      
      }
      
      
    }
    
  } else  if (y == "prec"){
    
    if (cpp){
      Jac[meanPart,meanPart] <- d_sigma_kappa_cpp(L = Lmat, D = dots$D, sigma = dots$sigma_y)
    } else {
      Jac[meanPart,meanPart] <- d_sigma_kappa(L = Lmat, D = dots$D, sigma = dots$sigma_y)      
    }
    
    
  } else if (y == "cor"){
    # Corelation matrix:
    corPart <- seq_len(nvar*(nvar-1)/2)
    
    if (cpp){
      Jac[meanPart,corPart] <- d_sigma_rho_cpp(L = Lmat, SD = dots$SD_y, A = dots$A, Dstar = dots$Dstar)
    } else {
      Jac[meanPart,corPart] <- d_sigma_rho(L = Lmat, SD = dots$SD_y, A = dots$A, Dstar = dots$Dstar)      
    }
    
    
    if (!metacor){
      sdPart <- max(corPart) + seq_len(nvar)  
      
      if (cpp){
        Jac[meanPart,sdPart] <- d_sigma_SD_cpp(L = Lmat, SD_IplusRho = dots$SD_IplusRho_y, In = dots$In, A = dots$A)
      } else {
        Jac[meanPart,sdPart] <- d_sigma_SD(L = Lmat, SD_IplusRho = dots$SD_IplusRho_y, In = dots$In, A = dots$A)
      }
      
    }
  }
  
  
  ### Add random effects (mostly same code):
  nEl <- nmod * (nmod+1) / 2
  if (randomEffects == "cov"){
    # Regular covs:
    Jac[varPart,varPart] <- as.matrix(Diagonal(nEl))
    
  } else if (randomEffects == "chol"){
    if (cpp){
      # Cholesky decomposition:
      Jac[varPart,varPart] <- d_sigma_cholesky_cpp(lowertri = dots$lowertri_randomEffects,L = dots$L_c, C = dots$C_c, 
                                                   In = dots$In_c)
    } else {
      # Cholesky decomposition:
      Jac[varPart,varPart] <- d_sigma_cholesky(lowertri = dots$lowertri_randomEffects,L = dots$L_c, C = dots$C_c, 
                                               In = dots$In_c)      
    }
    
  } else if (randomEffects == "ggm"){
    
    # Gaussian graphical model:
    netPart <- max(meanPart) + seq_len(nmod*(nmod-1)/2)
    scalingPart <- max(netPart) + seq_len(nmod)
    
    # if (metacor){
    #   
    #   if (cpp){
    #     Jac[varPart,netPart] <- d_sigma_omega_corinput_cpp(delta_IminOinv = dots$delta_IminOinv_randomEffects, 
    #                                                        L = dots$L_c, A = dots$A_c, delta = dots$delta_randomEffects,
    #                                                        Dstar = dots$Dstar_c, IminOinv = dots$IminOinv_randomEffects, 
    #                                                        In = dots$In_c)
    #   } else {
    #     Jac[varPart,netPart] <- d_sigma_omega_corinput(delta_IminOinv = dots$delta_IminOinv_randomEffects, 
    #                                                    L = dots$L_c, A = dots$A_c, delta = dots$delta_randomEffects,
    #                                                    Dstar = dots$Dstar_c, IminOinv = dots$IminOinv_randomEffects, 
    #                                                    In = dots$In_c)
    #   }
    #   
    # } else {
      if (cpp){
        Jac[varPart,netPart] <- d_sigma_omega_cpp(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                  A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c)
        # Jac[varPart,scalingPart] <- d_sigma_delta_cpp(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        #                                               In = dots$In_c, A = dots$A_c, delta = dots$delta_randomEffects)
        Jac[varPart,scalingPart] <- d_sigma_delta_cpp(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                      In = dots$In_c, A = dots$A_c)
        
      } else {
        Jac[varPart,netPart] <- d_sigma_omega(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                              A = dots$A_c, delta = dots$delta_randomEffects, Dstar = dots$Dstar_c)
        # Jac[varPart,scalingPart] <- d_sigma_delta(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
        #                                           In = dots$In_c, A = dots$A_c, delta = dots$delta_randomEffects)
        
        Jac[varPart,scalingPart] <- d_sigma_delta(L = dots$L_c, delta_IminOinv = dots$delta_IminOinv_randomEffects,
                                                  In = dots$In_c, A = dots$A_c)
        
      }
      
    # }
    
  } else  if (randomEffects == "prec"){
    
    if (cpp){
      Jac[varPart,varPart] <- d_sigma_kappa_cpp(L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects)
    } else {
      Jac[varPart,varPart] <- d_sigma_kappa(L = dots$L_c, D = dots$D_c, sigma = dots$sigma_randomEffects)
    }
    
  } else if (randomEffects == "cor"){
    # Corelation matrix:
    corPart <- max(meanPart) + seq_len(nmod*(nmod-1)/2)
    
    if (cpp){
      Jac[varPart,corPart] <- d_sigma_rho_cpp(L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c)
    } else {
      Jac[varPart,corPart] <- d_sigma_rho(L = dots$L_c, SD = dots$SD_randomEffects, A = dots$A_c, Dstar = dots$Dstar_c)      
    }

 
    # if (!metacor){
      sdPart <- max(corPart) + seq_len(nmod) 
      
      if (cpp){
        Jac[varPart,sdPart] <- d_sigma_SD_cpp(L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects, In = dots$In_c, 
                                          A = dots$A_c)
      } else {
        Jac[varPart,sdPart] <- d_sigma_SD(L = dots$L_c, SD_IplusRho = dots$SD_IplusRho_randomEffects, In = dots$In_c, 
                                          A = dots$A_c)        
      }

    # }
  }
  
  ####
  
  # Make sparse if needed:
  # Jac <- as(Jac, "Matrix")
  Jac <- sparseordense(Jac)
  
  # Return jacobian:
  return(Jac)
}

# Now for the full group:
d_phi_theta_meta_varcov <- function(prep){
  # model is already prepared!
  
  # d_phi_theta per group:
  d_per_group <- lapply(prep$groupModels,do.call,what=d_phi_theta_meta_varcov_group)
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",d_per_group)
}


