prepareModel <- function(x, model){
  # What model?
  framework <- model@model
  
  # Get the function:
  prepFun <- switch(framework,
        "lnm" = prepare_lnm,
        "ggm" = prepare_ggm,
        "rnm" = prepare_rnm,
        "gvar" = prepare_gvar,
        "varcov" = prepare_varcov,
        "cholesky" = prepare_cholesky)
# prepare:
  prep <- prepFun(x, model)
    
  # If the estimator is FIML, add the raw data:
  if (model@estimator == "FIML"){
    # Add the raw data to each group:
    for (g in seq_along(prep$groupModels)){
      prep$groupModels[[g]]$data <- model@sample@data[[g]]
    }
  }
  
 # Return:
  return(prep)
}