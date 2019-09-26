# Implied model for precision. Requires appropriate model matrices:
implied_meta_varcov <- function(model, all = FALSE){
  x <- formModelMatrices(model)

  # Implied covariance structures:
  x <- impliedcovstructures(x,type = model@types$y, all = all, name = "y")
  x <- impliedcovstructures(x,type = model@types$randomEffects, name = "randomEffects", all = all)

  for (g in seq_along(x)){
    
    # the 'meanstructure' is the varcov structure:
    sigma_y <- x[[g]]$sigma_y
    if (model@sample@corinput){
      x[[g]]$mu <- Vech(cov2cor(sigma_y), FALSE)
    } else {
      x[[g]]$mu <- Vech(sigma_y, TRUE)
    }
    
    # Form the var-cov matrix:
    x[[g]]$sigma <- x[[g]]$sigma_randomEffects + model@extramatrices$V
    x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
  }
 
  
  x
}
