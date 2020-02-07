impliedModel <- function(model, types, all = FALSE){
  # What model?
  framework <- model@model

  # Get the function:
  impFun <- switch(framework,
        "lvm" = implied_lvm,
        "varcov" = implied_varcov,
        # "gvar" = implied_gvar_norawts,
        "var1" = implied_var1,
        # "panelvar1" = implied_panelvar1,
        "dlvm1" = implied_dlvm1,
        "tsdlvm1" = implied_tsdlvm1,
        "meta_varcov" = implied_meta_varcov,
        "Ising" = implied_Ising,
        "ml_lvm" = implied_ml_lvm
        )
  
# implied:
  # mats <- formModelMatrices(model)
  imp <- impFun(model, all)
    
 # Return:
  return(imp)
}