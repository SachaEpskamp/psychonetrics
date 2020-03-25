prepareModel <- function(x, model){
  # What model?
  framework <- model@model
  
  # Get the function:
  if (model@cpp){
    prepFun <- switch(framework,
                      # "lnm" = prepare_lnm,
                      # "ggm" = prepare_ggm,
                      # "rnm" = prepare_rnm,
                      # "gvar" = prepare_gvar,
                      "varcov" = prepare_varcov_cpp, # <- Updated!
                      "lvm" = prepare_lvm_cpp, # <- updated!
                      "var1" = prepare_var1_cpp, # <- Updated!
                      # "panelvar1" = prepare_panelvar1,
                      "dlvm1" = prepare_dlvm1_cpp, # <- Updayed!
                      "tsdlvm1" = prepare_tsdlvm1_cpp, # <- Updated!
                      "meta_varcov" = prepare_meta_varcov_cpp, # <- updated!
                      "Ising" = prepare_Ising_cpp, # <- Updated!
                      "ml_lvm" = prepare_ml_lvm_cpp # <- Udated!!!
                      # "cholesky" = prepare_cholesky
    )    
  } else {
    prepFun <- switch(framework,
                      # "lnm" = prepare_lnm,
                      # "ggm" = prepare_ggm,
                      # "rnm" = prepare_rnm,
                      # "gvar" = prepare_gvar,
                      "varcov" = prepare_varcov,
                      "lvm" = prepare_lvm,
                      "var1" = prepare_var1,
                      # "panelvar1" = prepare_panelvar1,
                      "dlvm1" = prepare_dlvm1,
                      "tsdlvm1" = prepare_tsdlvm1,
                      "meta_varcov" = prepare_meta_varcov,
                      "Ising" = prepare_Ising,
                      "ml_lvm" = prepare_ml_lvm
                      # "cholesky" = prepare_cholesky
    )
  }

# prepare:
  prep <- prepFun(x, model)
    
  # If the estimator is FIML, add the raw data:
  if (model@estimator == "FIML"){
    # Add the raw data to each group:
    for (g in seq_along(prep$groupModels)){
      prep$groupModels[[g]]$fimldata <- model@sample@fimldata[[g]]
      prep$groupModels[[g]]$fulln <- model@sample@groups$nobs[[g]]
    }
  }
  
  # FIXME: Add the estimator to group modes for DWLS:
  for (g in seq_along(prep$groupModels)){
    prep$groupModels[[g]]$estimator <- model@estimator
  }
  
  # FIXME: Add Cpp to prep model:
  prep$cpp <- model@cpp
  
  # FIXME: Add fullFIML to the prep model:
  prep$fullFIML <- model@sample@fullFIML
  
  # FIXME Add WLS.W:
  if (model@estimator %in%  c("WLS","DWLS","ULS")){
    for (g in seq_along(prep$groupModels)){
      prep$groupModels[[g]]$WLS.W <- model@sample@WLS.W[[g]]
    }
  }
  
  # Add number of parameters:
  prep$nParFull = nrow(model@parameters)
  prep$nParFree = max(model@parameters$par[!is.na(model@parameters$par)])
  
 # Return:
  return(prep)
}