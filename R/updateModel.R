updateModel <- function(x, model){
  # New model:
  newMod <- model
  
  # Overwrite the current estimates:
  for (i in seq_along(x)){
    newMod@parameters$est[newMod@parameters$par==i] <- x[i]
  }
  newMod@modelmatrices <- formModelMatrices(newMod)
  newMod
}
