prepareModel <- function(x, model){
  # What model?
  framework <- model@model
  
  # Get the function:
  prepFun <- switch(framework,
        # "lnm" = prepare_lnm,
        # "ggm" = prepare_ggm,
        # "rnm" = prepare_rnm,
        "gvar" = prepare_gvar,
        "varcov" = prepare_varcov,
        "lvm" = prepare_lvm,
        "var1" = prepare_var1
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
  
 # Return:
  return(prep)
}