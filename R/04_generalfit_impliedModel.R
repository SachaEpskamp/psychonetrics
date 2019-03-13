impliedModel <- function(model, types, all = FALSE){
  # What model?
  framework <- model@model

  # Get the function:
  impFun <- switch(framework,
        "lvm" = implied_lvm,
        "varcov" = implied_varcov,
        "gvar" = implied_gvar_norawts)
  
# implied:
  # mats <- formModelMatrices(model)
  imp <- impFun(model, all)
    
 # Return:
  return(imp)
}