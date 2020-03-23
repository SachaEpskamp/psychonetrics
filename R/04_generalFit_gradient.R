# General gradient function:
psychonetrics_gradient <- function(x, model){
  # I need an estimator part, a model part and a manual part
  # Prepare
  # message("Prep model...")
  # Prepare model:
  
  # if (model@cpp){
  #   prep <- prepareModel_cpp(x, model) # <- upated!
  # } else {
    prep <- prepareModel(x, model)  
  # }
  


  # estimator part:
  if (model@cpp){
    
    estimatorJacobian <- switch(
      model@estimator,
      "ML" = switch(model@distribution,
                    "Gaussian" = jacobian_gaussian_sigma_cpp, # <- Updated!
                    "Ising" = jacobian_Ising_cpp # <- updated!
      ),
      "ULS" = switch(model@distribution,
                     "Gaussian" = ULS_gradient_Gauss_cpp),  # <- updated!
      "WLS" = switch(model@distribution,
                     "Gaussian" = ULS_gradient_Gauss_cpp),  # <- updated!
      "DWLS" = switch(model@distribution,
                      "Gaussian" = ULS_gradient_Gauss_cpp), # <- updated!
      "FIML" = switch(model@distribution,
                      "Gaussian" = jacobian_fiml_gaussian_sigma_cpp) # <- Updated!
    )
    
    
  } else {
    
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
    
  }
 

  # message("Estimator part...")
  estimatorPart <- estimatorJacobian(prep)
  
  # model part:
  if (model@cpp){
    modelJacobian <- switch(
      model@model,
      # "lnm" = d_phi_theta_lnm,
      # "ggm" = d_phi_theta_ggm,
      # "rnm" = d_phi_theta_rnm,
      # "gvar" = ifelse(model@rawts,d_phi_theta_gvar_rawts,d_phi_theta_gvar),
      "varcov" = d_phi_theta_varcov_cpp, # <- updated!
      "lvm" = d_phi_theta_lvm_cpp, # <- updated!
      "var1" = d_phi_theta_var1_cpp, # <- updated!
      # "panelvar1" = d_phi_theta_panelvar1,
      "dlvm1" = d_phi_theta_dlvm1_cpp, # <- updated!
      "tsdlvm1" = d_phi_theta_tsdlvm1_cpp, # <- updated!
      "meta_varcov" = d_phi_theta_meta_varcov_cpp, # <- updated!
      "Ising" = d_phi_theta_Ising_cpp, # <- updated!
      "ml_lvm" = d_phi_theta_ml_lvm_cpp # <- updated!
      # "cholesky" = d_phi_theta_cholesky
    )    
  } else {
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
  }

  

  # message("Model part...")
  modelPart <- sparseordense(modelJacobian(prep))
  
  

  # Manual part:
  if (model@cpp){
    manualPart <- Mmatrix_cpp(model@parameters)
  } else {
    manualPart <- Mmatrix(model@parameters)  
  }
 

  # FIXME: partial cpp
  if (model@cpp){
    if (is(modelPart, "sparseMatrix")){
      Jac <- gradient_inner_cpp_DSS(as.matrix(estimatorPart), modelPart, manualPart)
    } else {
      Jac <- gradient_inner_cpp_DDS(as.matrix(estimatorPart), as.matrix(modelPart), manualPart)
    }
    
  } else {

    # Full Jacobian:
    Jac <- estimatorPart %*% modelPart %*% manualPart    
  }

  
  # Return:
  return(as.vector(Jac))
}