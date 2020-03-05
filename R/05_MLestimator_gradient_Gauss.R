jacobian_gaussian_group_sigmaVersion_meanPart <- function(sigma,mu,means,kappa,...){
  # Mean part:
  # grad_mean <- -2 * t(means - mu) %*% kappa
  grad_mean <- -2 * t(means - mu) %*% kappa
  grad_mean
}

jacobian_gaussian_group_sigmaVersion_sigmaPart <- function(S,means,mu,sigma,D,kappa,...){
    # if (nrow(Drawts) != ncol(Drawts)){
    #   D <- Diagonal(n = nrow(sigma)^2)
    # } 
    # sigma part:
    # grad_sigma <- - t(Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma)) %*% (kappa %x% kappa) %*% D
  
  # Try fancier method to avoid memory problems...
  n <- ncol(S)
  # mat <- Matrix(as.vector(Vec(S) + Vec((means - mu) %*% t(means - mu)) - Vec(sigma)),n,n)
  mat <- S + (means - mu) %*% t(means - mu) - sigma
  grad_sigma <- t(-t(D) %*% Vec(kappa %*% mat %*% kappa))

  as.matrix(grad_sigma)
}


# jacobian function per group:
jacobian_gaussian_group_sigma <- function(...,Drawts,mu,sigma, meanstructure = TRUE, corinput = FALSE){

  
  # sigma part:
  grad_sigma <- jacobian_gaussian_group_sigmaVersion_sigmaPart(mu=mu,sigma=sigma,...,Drawts=Drawts)

  # FIXME: Nicer to not bother with computing the parts not needed
  # Cut out variances if needed:
  if (corinput){
    keep <- diag(ncol(sigma))[lower.tri(diag(ncol(sigma)),diag=TRUE)] != 1
    grad_sigma <- as(grad_sigma[,keep, drop=FALSE], "matrix")
  }
  
  # Cut out means if needed:
  if (meanstructure){
    # Mean part:
    grad_mean <- jacobian_gaussian_group_sigmaVersion_meanPart(mu=mu,sigma=sigma,...)
    
    # Put int output:
    Out <- cbind(grad_mean,grad_sigma)
  } else {
    Out <- grad_sigma
  }
  
  # if (missing(Drawts)){
  #   Drawts <- Diagonal(ncol(Out))
  # }
  
  # Combine and return:
  Out# %*% Drawts
}

# Now for all groups:
jacobian_gaussian_sigma <- function(prep){
  # model is already prepared!

  # d_phi_theta per group:
  g_per_group <- lapply(prep$groupModels,do.call,what=jacobian_gaussian_group_sigma)
  
  # Weight:
  for (i in 1:length(prep$groupModels)){
    g_per_group[[i]] <- (prep$nPerGroup[i] / prep$nTotal) * g_per_group[[i]]
  }
  
  # FIXME: Computationall it is nicer to make the whole object first then fill it
  # Bind by colum and return: 
  Reduce("cbind",g_per_group)
}
