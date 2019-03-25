# Latent network model creator
ggm <- function(
  ...
){
 model <- varcov(...,type = "ggm")
  
  # Return model:
  return(model)
}
