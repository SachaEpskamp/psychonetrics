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
  }
  
    
  newMod
}
