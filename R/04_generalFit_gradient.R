# General gradient function:
psychonetrics_gradient <- function(x, model){
  # I need an estimator part, a model part and a manual part
  
  # Prepare
  prep <- prepareModel(x, model)

  # estimator part:
  estimatorJacobian <- switch(
    model@estimator,
    "ML" = switch(model@distribution,
                  "Gaussian" = jacobian_gaussian_sigma
    ),
    "ULS" = switch(model@distribution,
                "Gaussian" = ULS_gradient_Gauss)
  )
  estimatorPart <- estimatorJacobian(prep)
  
  # model part:
  modelJacobian <- switch(
    model@model,
    "lnm" = d_phi_theta_lnm,
    "ggm" = d_phi_theta_ggm,
    "rnm" = d_phi_theta_rnm
  )
  modelPart <- modelJacobian(prep)
 
  # Manual part:
  manualPart <- Mmatrix(model@parameters)
  
  # Full Jacobian:
  Jac <- estimatorPart %*% modelPart %*% manualPart

  # Return:
  return(as.vector(Jac))
}