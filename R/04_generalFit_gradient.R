# General gradient function:
psychonetrics_gradient <- function(x, model){
  # I need an estimator part, a model part and a manual part
  # Prepare
  # message("Prep model...")
  prep <- prepareModel(x, model)

  # estimator part:
  estimatorJacobian <- switch(
    model@estimator,
    "ML" = switch(model@distribution,
                  "Gaussian" = jacobian_gaussian_sigma,
                  "Ising" = jacobian_Ising
    ),
    "ULS" = switch(model@distribution,
                "Gaussian" = ULS_gradient_Gauss),
    "WLS" = switch(model@distribution,
                   "Gaussian" = ULS_gradient_Gauss),
    "DWLS" = switch(model@distribution,
                   "Gaussian" = ULS_gradient_Gauss),
    "FIML" = switch(model@distribution,
                   "Gaussian" = jacobian_fiml_gaussian_sigma)
  )

  # message("Estimator part...")
  estimatorPart <- estimatorJacobian(prep)
  
  # model part:
  modelJacobian <- switch(
    model@model,
    # "lnm" = d_phi_theta_lnm,
    # "ggm" = d_phi_theta_ggm,
    # "rnm" = d_phi_theta_rnm,
    # "gvar" = ifelse(model@rawts,d_phi_theta_gvar_rawts,d_phi_theta_gvar),
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

  # message("Model part...")
  modelPart <- modelJacobian(prep)
 
  # Manual part:
  manualPart <- Mmatrix(model@parameters)

  # Full Jacobian:
  Jac <- estimatorPart %*% modelPart %*% manualPart

  # Return:
  return(as.vector(Jac))
}