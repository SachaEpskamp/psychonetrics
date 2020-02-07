prepareModel <- function(x, model){
  # What model?
  framework <- model@model
  
  # Get the function:
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
  
  # FIXME Add WLS.V:
  if (model@estimator %in%  c("WLS","DWLS","ULS")){
    for (g in seq_along(prep$groupModels)){
      prep$groupModels[[g]]$WLS.V <- model@sample@WLS.V[[g]]
    }
  }
  
 # Return:
  return(prep)
}