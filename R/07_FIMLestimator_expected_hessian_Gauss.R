# expected_hessian_fiml_Gaussian_group_meanPart <- function(kappa,data,sigma,...){
#   n <- nrow(data)
#   
#   # Spectral shift kappa (will do nothing if kappa is pos def):
#   kappa <- spectralshift(kappa)
#   
#   # Current jacobian:
#   curFish <- Matrix(0,nrow=nrow(kappa),ncol=nrow(kappa))
# 
#   # For every subject:
#   for (i in seq_len(n)){
#     if (!any(is.na(data[i,]))){
#       kappa_p <- kappa
#       # Elimination matrix:
#       L <- Diagonal(n=nrow(kappa))
#     } else {
#       obs <- unlist(!is.na(data[i,]))
#       
#       # Skip if nothing to do:
#       if (!any(obs)) {
#         next
#       }
#       
#       sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
#       kappa_p <- corpcor::pseudoinverse(sig_p)
#       
#       # Handle possible non positive definiteness:
#       kappa_p <- spectralshift(kappa_p)
#       
#       # Find the proper elimination matrix:
#       inds <- which(obs)
#       L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(sigma)))
#     }
#     
#     curFish <- curFish + t(L) %*% kappa_p %*% L
#     
#   }
#   
#   # Mean part:
#   (2/n) * curFish
# }
# 
# expected_hessian_fiml_Gaussian_group_varPart <- function(kappa,D,data,sigma,...){
#   
#   n <- nrow(data)
#   
#   # Spectral shift kappa (will do nothing if kappa is pos def):
#   kappa <- spectralshift(kappa)
#   
#   # DUmmy sigma for indices:
#   dumSig <- matrix(1:(ncol(sigma)^2),ncol(sigma),ncol(sigma))
#   
#   # Current expHEssian:
#   curFish <- Matrix(0,nrow=nrow(sigma) * (nrow(sigma)+1) / 2,ncol=nrow(sigma) * (nrow(sigma)+1) / 2)
#   
#   # Kappa kronecker kappa for full cases:
#   KkronK <- (kappa %(x)% kappa)
#   
#   # For every subject:
#   for (i in seq_len(n)){
#     
#     if (!any(is.na(data[i,]))){
#       sigma_p <- sigma
#       # Elimination matrix:
#       L <- Diagonal(n=nrow(sigma)^2)
#       KkronK_p <- KkronK
#     } else {
#       obs <- unlist(!is.na(data[i,]))
#       
#       # Skip if nothing to do:
#       if (!any(obs)) {
#         next
#       }
#       
#       sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]
#       kappa_p <- corpcor::pseudoinverse(sig_p)
#       
#       # Handle possible non positive definiteness:
#       kappa_p <- spectralshift(kappa_p)
#       
#       # Kronecker:
#       KkronK_p <-  (kappa_p %(x)% kappa_p)
#       
#       # Find the proper elimination matrix:
#       inds <- c(dumSig[obs,obs,drop=FALSE])
#       L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(sigma)^2))
#     }
#     
#     # Update Fisher:
#     curFish <- curFish + t(D) %*% t(L) %*% KkronK_p %*% L %*% D
#   }
#   
#   # Return:
#   (1/n) *curFish
# }

# Per group:
expected_hessian_fiml_Gaussian_subgroup <- function(dat,sigma,kappa,mu,...){
  obs <- !as.vector(dat$pattern)
  
  sig_p <- as.matrix(sigma)[obs,obs,drop=FALSE]

  kappa_p <- corpcor::pseudoinverse(sig_p)
  
  # Handle possible non positive definiteness:
  kappa_p <- spectralshift(kappa_p)
  
  # Log determinant:
  # logdetK_p <- log(det(kappa_p))
  
  # Means:
  mu_p <- mu[obs,]
  dat$n * t(dat$L) %*% expected_hessian_Gaussian_group(S = dat$S,kappa = kappa_p, mu = mu_p, means = dat$means,sigma = sig_p,D=dat$D) %*% dat$L
}

# Fit function per group:
expected_hessian_fiml_Gaussian_group <- function(fimldata,fulln,sigma,kappa,mu,...){
  
  # Subgroup models:
  1/fulln * Reduce("+", lapply(fimldata,expected_hessian_fiml_Gaussian_subgroup,sigma=sigma,kappa=kappa,mu=mu))
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
  exph_per_group <- lapply(prep$groupModels,do.call,what=expected_hessian_fiml_Gaussian_group)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    exph_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * exph_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("bdiag",exph_per_group)
}
