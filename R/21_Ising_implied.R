# Implied model for precision. Requires appropriate model matrices:
implied_Ising <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }
  
  if (model@types$beta_model == "log_beta"){
    # For each group:
    nGroup <- length(x)
    
    for (g in 1:nGroup){
      x[[g]]$beta <- exp(x[[g]]$log_beta)
    }
  }
 
  
  x
}
