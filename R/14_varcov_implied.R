# Implied model for precision. Requires appropriate model matrices:
implied_varcov <- function(model, all = FALSE){
  x <- formModelMatrices(model)
  
  # Implied covariance structures:
  x <- impliedcovstructures(x,type = model@types$y, all = all)

  for (g in seq_along(x)){
    if (is.null(x[[g]]$kappa)){
      x[[g]]$kappa <- solve_symmetric(x[[g]]$sigma, logdet = TRUE)
    } else {
      if (is.null(attr(x[[g]]$kappa, "logdet"))){
        ev <- eigen(x[[g]]$kappa, symmetric=TRUE, only.values=TRUE)$values
        
        # Check pos-def:
        if(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0) {
            attr(x[[g]]$kappa, "logdet") <- log(.Machine$double.eps)  
          
        } else {
          attr(x[[g]]$kappa, "logdet") <- log(det(x[[g]]$kappa))
        }
      }
    }
  }
 
  
  x
}
