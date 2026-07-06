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
                    "Spin" = jacobian_Ising_cpp, # <- updated!
                    "TwoLevelGaussian" = jacobian_gaussian2L_sigma_cpp # <- updated!
      ),
      "PML" = switch(model@distribution,
                     "Gaussian" = jacobian_gaussian_sigma_cpp,
                     "Spin" = jacobian_Ising_cpp
      ),
      "PFIML" = switch(model@distribution,
                       "Gaussian" = jacobian_fiml_gaussian_sigma_cpp
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
                    "Spin" = jacobian_Ising,
                    "TwoLevelGaussian" = jacobian_gaussian2L_sigma
      ),
      "PML" = switch(model@distribution,
                     "Gaussian" = jacobian_gaussian_sigma,
                     "Spin" = jacobian_Ising
      ),
      "PFIML" = switch(model@distribution,
                       "Gaussian" = jacobian_fiml_gaussian_sigma
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
      "panelvar" = d_phi_theta_panelvar_cpp,
      "tsdlvm1" = d_phi_theta_tsdlvm1_cpp, # <- updated!
      "meta_varcov" = d_phi_theta_meta_varcov_cpp, # <- updated!
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising_cpp, # <- updated!
      "ml_lvm" = d_phi_theta_ml_lvm_cpp, # <- updated!
      "meta_lvm" = d_phi_theta_meta_lvm_cpp,
      "meta_var1" = d_phi_theta_meta_var1_cpp
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
      "panelvar" = d_phi_theta_panelvar,
      "ml_var1" = d_phi_theta_ml_var1,
      "tsdlvm1" = d_phi_theta_tsdlvm1,
      "meta_varcov" = d_phi_theta_meta_varcov,
      "Ising" = ,
      "BlumeCapel" = d_phi_theta_Ising,
      "ml_lvm" = d_phi_theta_ml_lvm,
      "meta_lvm" = d_phi_theta_meta_lvm,
      "meta_var1" = d_phi_theta_meta_var1
      # "cholesky" = d_phi_theta_cholesky
    )
  }



  # message("Model part...")
  # Model part. For the gradient the estimator part is a single row vector,
  # so the dense matrix products below are already fast; converting a merely
  # sparse (but not diagonal) model Jacobian to a sparse Matrix costs more
  # than the sparse products save. Hence, only exactly diagonal Jacobians
  # (cheap to store and multiply) take the sparse path here. The Fisher
  # information (04_generalfit_FisherInformation.R) does use the full
  # sparse-or-dense classification, where the sparse path is a large win:
  modelPart <- modelJacobian(prep)
  if (!is.matrix(modelPart)) modelPart <- as.matrix(modelPart)
  if (diag_sparse_dense_cpp(modelPart) == 0L){
    modelPart <- as(modelPart, "dMatrix")
  }
  
  

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

  # Compute gradient:
  grad <- as.vector(Jac)

  # Add penalty gradient for PML/PFIML estimator:
  if (model@estimator %in% c("PML", "PFIML")) {
    grad <- addPenaltyGradient(grad, x, model)
  }

  # Return:
  return(grad)
}