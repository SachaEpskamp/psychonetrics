numeric_FisherInformation <- function(model){
  model <- expectedmodel(model)
  # 2 * sum(model@sample@groups$nobs) * numDeriv::jacobian(psychonetrics_gradient,parVector(model),model=model)
  
  # Unit information instead:
  0.5 * numDeriv::jacobian(psychonetrics_gradient,parVector(model),model=model)
}

# General Fisher information former:
psychonetrics_FisherInformation <- function(model, analytic = TRUE){
  # Maybe do numeric?
  if (!analytic){
    return(numeric_FisherInformation(model))
  }
  
  # Prepare (FIXME: x not needed but I am laxy...)
  prep <- prepareModel(parVector(model), model)
  
  # # estimator part Jacobian:
  # estimatorJacobian <- switch(
  #   model@estimator,
  #   "ML" = switch(model@distribution,
  #                 "Gaussian" = jacobian_gaussian_sigma
  #   )
  # )
  # estimatorPartJacobian <- estimatorJacobian(prep)
  
  # estimator part expected Hessian:
  estimatorHessian <- switch(
    model@estimator,
    "ML" = switch(model@distribution,
                  "Gaussian" = expected_hessian_Gaussian,
                  "Ising" = expected_hessian_Ising
    ),
    "ULS" = switch(model@distribution,
                  "Gaussian" = expected_hessian_ULS_Gaussian
    ),
    "WLS" = switch(model@distribution,
                   "Gaussian" = expected_hessian_ULS_Gaussian
    ),
    "DWLS" = switch(model@distribution,
                   "Gaussian" = expected_hessian_ULS_Gaussian
    ),
    "FIML" = switch(model@distribution,
                   "Gaussian" = expected_hessian_fiml_Gaussian
    )
  )
  estimatorPartHessian <- estimatorHessian(prep)
  
  # model part:
  modelJacobian <- switch(
    model@model,
    # "lnm" = d_phi_theta_lnm,
    # "ggm" = d_phi_theta_ggm,
    # "rnm" = d_phi_theta_rnm,
    # "gvar" =  ifelse(model@rawts,d_phi_theta_gvar_rawts,d_phi_theta_gvar),
    "varcov" = d_phi_theta_varcov,
    "lvm" = d_phi_theta_lvm,
    "var1" = d_phi_theta_var1,
    # "panelvar1" = d_phi_theta_panelvar1,
    "dlvm1" = d_phi_theta_dlvm1,
    "tsdlvm1" = d_phi_theta_tsdlvm1,
    "meta_varcov" = d_phi_theta_meta_varcov,
    "Ising" = d_phi_theta_Ising,
    "ml_lvm" = d_phi_theta_ml_lvm
    # "cholesky" = d_phi_theta_cholesky
  )
  modelPart <- modelJacobian(prep)
  
  # Manual part:
  manualPart <- Mmatrix(model@parameters)

  # Compute fisher information and return:
  # Fisher <- 2 * prep$nTotal * t(manualPart) %*% t(modelPart) %*% estimatorPartHessian %*% modelPart %*% manualPart

  # Unit information instead:
  Fisher <- 0.5 * t(manualPart) %*% t(modelPart) %*% estimatorPartHessian %*% modelPart %*% manualPart
  
  as.matrix(Fisher)
}
# 
# 
# # Fisher information for LNM model
# Fisher_lnm <- function(model){
#   # Prepare
#   prep <- prepare_lnm(parVector(model), model)
# 
#   # d_phi_theta:
#   d_phi_theta <- d_phi_theta_lnm(prep)
#   
#   # Model matrix:
#   M <- Mmatrix(model@parameters)
#   
#   # Number of groups:
#   nGroups <- nrow(model@sample@groups)
#   
#   # Total number of parameters:
#   nTotal <- nrow(M)
#   nPar_perGroup <- nTotal / nGroups
#   
#   # number of constrained parameters:
#   nConstrained <- ncol(M)
#   
#   # Number of variables:
#   nVars <- nrow(model@sample@variables)
#   
#   # n means:
#   nMeans <- nVars
#   
#   # N var/covs:
#   nVarCov <- nVars*(nVars+1)/2
#   
#   # total per group:
#   nPerGroup <- nMeans + nVarCov
#   
#   # Create empty fisher information matrix:
#   I <- Matrix(0,nTotal,nTotal)
#   # For each group
#   for (g in 1:nGroups){
# 
#     # Indices of the mean part:
#     meanPart <- (g-1)*nPerGroup + seq_len(nMeans)
#     varPart <- (g-1)*nPerGroup + nMeans + seq_len(nVarCov)
#     groupPart <- (g-1)*nPar_perGroup + seq_len(nPar_perGroup)
#     kappa <- prep$groupModels[[g]]$kappa
#     
#     # Duplication mat:
#     D <- prep$groupModels[[g]]$D
# 
#     # Fill in
#     I[(g-1)*nPar_perGroup + seq_len(nPar_perGroup),(g-1)*nPar_perGroup + seq_len(nPar_perGroup)] <-
#       # (prep$nPerGroup[g] / prep$nTotal) *
#       ( prep$nTotal / prep$nPerGroup[g]) *
#       2*t(d_phi_theta[meanPart,groupPart]) %*% kappa %*% d_phi_theta[meanPart,groupPart] +
#       t(d_phi_theta[varPart,groupPart]) %*% t(D) %*% (kappa %x% kappa) %*% D %*% d_phi_theta[varPart,groupPart]
# 
# # 
# #     
# #     # Fill the relevant elements:
# #     # for (i in (g-1)*nPar_perGroup + seq_len(nPar_perGroup)){
# #     #   # Derivatives of sigma for i:
# #     #   dSigmai <- sparseMatrix(
# #     #     i = row(kappa)[lower.tri(kappa,diag=TRUE)],
# #     #     j = col(kappa)[lower.tri(kappa,diag=TRUE)],
# #     #     x = d_phi_theta[varPart, i], symmetric = TRUE)
# #     #   
# #     #   for (j in ((g-1)*nPar_perGroup + 1):i){
# #     #     # Derivatives of sigma for j:
# #     #     dSigmaj <- sparseMatrix(
# #     #       i = row(kappa)[lower.tri(kappa,diag=TRUE)],
# #     #       j = col(kappa)[lower.tri(kappa,diag=TRUE)],
# #     #       x = d_phi_theta[varPart, j], symmetric = TRUE)
# #     #     
# #     #     I[i,j] <- I[j,i] <- 2 * t(d_phi_theta[meanPart,i,drop=FALSE]) %*% kappa %*% d_phi_theta[meanPart,j,drop=FALSE] + 
# #     #        sum(diag(kappa %*% dSigmai %*% kappa %*% dSigmaj ))
# #     #   }
# #     # }
#   }
#   
# 
#   
#   # Return fisher information:
#   Fish <- t(M) %*% I %*% M
#   
#   # Return:
#   return(Fish)
# }