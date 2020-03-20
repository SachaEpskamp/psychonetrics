# Implied model for precision. Requires appropriate model matrices:
implied_Ising <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }
  
  x
}
