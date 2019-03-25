# Latent network model creator
cholesky <- function(
  ...
){
  model <- varcov(...,type = "chol")
  
  # Return model:
  return(model)
}
