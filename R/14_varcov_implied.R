# Implied model for precision. Requires appropriate model matrices:
implied_varcov <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }

  # Implied covariance structures:
  if (model@cpp){
    x <- impliedcovstructures_cpp(x,type = model@types$y, all = all)  
  } else {
    x <- impliedcovstructures(x,type = model@types$y, all = all)
  }
  

  for (g in seq_along(x)){
    # Check mu (if it is not there, this is because meanstructure = FALSE):
    if (is.null(x[[g]]$mu)){
      x[[g]]$mu <- as.matrix(model@sample@means[[g]])
    }

    # Check kappa
    if (is.null(x[[g]]$kappa)){
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    } else {
      if (is.null(attr(x[[g]]$kappa, "logdet"))){
        # ev <- eigen(x[[g]]$kappa, symmetric=TRUE, only.values=TRUE)$values
        
        # Check pos-def:
        # if(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0) {
        if (!sympd_cpp(x[[g]]$kappa)){
            attr(x[[g]]$kappa, "logdet") <- log(.Machine$double.eps)  
        } else {
          attr(x[[g]]$kappa, "logdet") <- determinant(x[[g]]$kappa, logarithm = TRUE)$modulus
        }
      }
    
    }
    

  }
 
  
  x
}
