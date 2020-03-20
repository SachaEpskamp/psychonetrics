# Implied model for precision. Requires appropriate model matrices:
implied_meta_varcov <- function(model, all = FALSE){

  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }

  
  if (model@cpp){
    
    # Implied covariance structures:
    x <- impliedcovstructures_cpp(x,type = model@types$y, all = all, name = "y")
    x <- impliedcovstructures_cpp(x,type = model@types$randomEffects, name = "randomEffects", all = all)
    
  } else {
    # Implied covariance structures:
    x <- impliedcovstructures(x,type = model@types$y, all = all, name = "y")
    x <- impliedcovstructures(x,type = model@types$randomEffects, name = "randomEffects", all = all)
  }
  


  for (g in seq_along(x)){
    
    est <- model@extramatrices$Vestimation
   
    # est <- model@extramatrices$Vmethod
    
    if (est == "averaged"){
    # if (est == "pooled"){
      # the 'meanstructure' is the varcov structure:
      sigma_y <- x[[g]]$sigma_y
      if (model@sample@corinput){
        x[[g]]$mu <- Vech(cov2cor(sigma_y), FALSE)
      } else {
        x[[g]]$mu <- Vech(sigma_y, TRUE)
      }
      

      # Form the var-cov matrix:
      x[[g]]$sigma <- x[[g]]$sigma_randomEffects + model@extramatrices[['V']]
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    } else {
      nStudy <- model@sample@groups$nobs[g]
        
      # Per group estimation, do this per group:
      sigma_y <- x[[g]]$sigma_y
      if (model@sample@corinput){
        x[[g]]$mu <- lapply(seq_len(nStudy),function(x)Vech(cov2cor(sigma_y), FALSE))
      } else {
        x[[g]]$mu <- lapply(seq_len(nStudy),function(x)Vech(sigma_y, TRUE))
      }
      x[[g]]$mu <- lapply(x[[g]]$mu,as.vector)
      
      # Form the var-cov matrix:
      x[[g]]$sigma <- lapply(seq_len(nStudy),function(i) x[[g]]$sigma_randomEffects + model@extramatrices$Vall[[i]])
      x[[g]]$sigma <- lapply( x[[g]]$sigma, as.matrix)
      x[[g]]$kappa <- lapply(x[[g]]$sigma,solve_symmetric,logdet=TRUE)
      x[[g]]$kappa <- lapply( x[[g]]$kappa, as.matrix)
    }
    
    # # Kappa, sigma and mu never sparse:
    # x[[g]]$mu <- as.matrix(x[[g]]$mu)
    # x[[g]]$kappa <- as.matrix(x[[g]]$kappa)
    # x[[g]]$sigma <- as.matrix(x[[g]]$sigma)
  

  }
 
  
  x
}
