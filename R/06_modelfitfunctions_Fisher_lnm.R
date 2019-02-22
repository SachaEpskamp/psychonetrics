# Fisher information for LNM model
Fisher_lnm <- function(model){
  # Prepare
  prep <- prepare_lnm(parVector(model), model)

  # d_phi_theta:
  d_phi_theta <- d_phi_theta_lnm(prep)
  
  # Model matrix:
  M <- Mmatrix(model@parameters)
  
  # Number of groups:
  nGroups <- nrow(model@sample@groups)
  
  # Total number of parameters:
  nTotal <- nrow(M)
  nPar_perGroup <- nTotal / nGroups
  
  # number of constrained parameters:
  nConstrained <- ncol(M)
  
  # Number of variables:
  nVars <- nrow(model@sample@variables)
  
  # n means:
  nMeans <- nVars
  
  # N var/covs:
  nVarCov <- nVars*(nVars+1)/2
  
  # total per group:
  nPerGroup <- nMeans + nVarCov
  
  # Create empty fisher information matrix:
  I <- Matrix(0,nTotal,nTotal)
  # For each group
  for (g in 1:nGroups){
    # Indices of the mean part:
    meanPart <- (g-1)*nPerGroup + seq_len(nMeans)
    varPart <- (g-1)*nPerGroup + nMeans + seq_len(nVarCov)
    kappa <- prep$groupModels[[g]]$kappa
    
    D <- prep$groupModels[[g]]$D
    I[(g-1)*nPar_perGroup + seq_len(nPar_perGroup),(g-1)*nPar_perGroup + seq_len(nPar_perGroup)] <-
      2*t(d_phi_theta[meanPart,]) %*% kappa %*% d_phi_theta[meanPart,] +
      t(d_phi_theta[varPart,]) %*% t(D) %*% (kappa %(x)% kappa) %*% D %*% d_phi_theta[varPart,]

# 
#     
#     # Fill the relevant elements:
#     # for (i in (g-1)*nPar_perGroup + seq_len(nPar_perGroup)){
#     #   # Derivatives of sigma for i:
#     #   dSigmai <- sparseMatrix(
#     #     i = row(kappa)[lower.tri(kappa,diag=TRUE)],
#     #     j = col(kappa)[lower.tri(kappa,diag=TRUE)],
#     #     x = d_phi_theta[varPart, i], symmetric = TRUE)
#     #   
#     #   for (j in ((g-1)*nPar_perGroup + 1):i){
#     #     # Derivatives of sigma for j:
#     #     dSigmaj <- sparseMatrix(
#     #       i = row(kappa)[lower.tri(kappa,diag=TRUE)],
#     #       j = col(kappa)[lower.tri(kappa,diag=TRUE)],
#     #       x = d_phi_theta[varPart, j], symmetric = TRUE)
#     #     
#     #     I[i,j] <- I[j,i] <- 2 * t(d_phi_theta[meanPart,i,drop=FALSE]) %*% kappa %*% d_phi_theta[meanPart,j,drop=FALSE] + 
#     #        sum(diag(kappa %*% dSigmai %*% kappa %*% dSigmaj ))
#     #   }
#     # }
  }
  

  
  # Return fisher information:
  Fish <- t(M) %*% I %*% M
  
  # Return:
  return(Fish)
}