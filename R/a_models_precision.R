# Latent network model creator
precision <- prec <- function(
  ...
){
 model <- varcov(...,type = "prec")
  
  # Return model:
  return(model)
}
