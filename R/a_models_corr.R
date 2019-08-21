# Correlation functio:
corr <- function(
  ...
){
 model <- varcov(...,type = "cor")
  
  # Return model:
  return(model)
}
