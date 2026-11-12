updateModel <- function(x, model, updateMatrices = FALSE){
  # New model:
  newMod <- model
  # Overwrite the current estimates:
  for (i in seq_along(x)){
    newMod@parameters$est[newMod@parameters$par==i] <- x[i]
  }

  # Add implied too:
  if (updateMatrices){
    newMod@modelmatrices <- impliedModel(newMod, model@types, all = TRUE)  

    # Composite variable models: show the computed variances and intercepts of
    # the composites in the parameter table:
    if (newMod@model == "cvm"){
      newMod <- cvm_fill_computed(newMod)
    }
  }
  
    
  newMod
}
