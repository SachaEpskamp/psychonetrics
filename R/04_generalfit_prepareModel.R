prepareModel <- function(x, model){
  # What model?
  framework <- model@model
  
  # Get the function:
  prepFun <- switch(framework,
        "lnm" = prepare_lnm,
        "ggm" = prepare_ggm)
  
  # Run and return:
  prepFun(x, model)
}