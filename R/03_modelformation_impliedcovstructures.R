impliedcovstructures <- function(
  x, # List with matrices
 name = "",
  type = "cov",
 all = FALSE
){
  # For each group:
  nGroup <- length(x)
  
  if (name != ""){
    sigma = paste0("sigma","_",name)
    omega = paste0("omega","_",name)
    delta = paste0("delta","_",name)
    kappa = paste0("kappa","_",name)
    lowertri = paste0("lowertri","_",name)
    IminOinv = paste0("IminOinv","_",name)
    delta_IminOinv = paste0("delta_IminOinv","_",name)
  } else {
    sigma = "sigma"
    omega = "omega"
    delta = "delta"
    kappa = "kappa"
    lowertri = "lowertri"
    IminOinv = "IminOinv"
    delta_IminOinv = "delta_IminOinv"
  }
  
  for (g in seq_len(nGroup)){
    # Form the models:
    # Contemporaneous:
    if (type == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- as(solve_symmetric(x[[g]][[sigma]]), "Matrix")
          x[[g]][[omega]]  <- as(qgraph::wi2net(as.matrix(x[[g]][[kappa]])),"Matrix")          
        }
      }
    } else if(type == "chol"){
      # form cov matrix:
      x[[g]][[sigma]] <- as(x[[g]][[lowertri]] %*% t(x[[g]][[lowertri]]), "Matrix")
      
      # Return precision and network if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- as(solve_symmetric(x[[g]][[sigma]]), "Matrix")
          x[[g]][[omega]]  <- as(qgraph::wi2net(as.matrix(x[[g]][[kappa]])),"Matrix")
        }

      }
    } else if (type == "ggm"){
      x[[g]][[sigma]] <- x[[g]][[delta]] %*% solve_symmetric(Diagonal(ncol(x[[g]][[omega]])) - x[[g]][[omega]]) %*% x[[g]][[delta]]
      
      # Stuff needed if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- solve_symmetric(x[[g]][[sigma]])
        }
        
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]][[IminOinv]] <- as(solve_symmetric(Diagonal(ncol(x[[g]][[omega]])) - x[[g]][[omega]]), "Matrix")
        x[[g]][[delta_IminOinv]] <- x[[g]][[delta]] %*% x[[g]][[IminOinv]]
      }
    } else if (type == "prec"){
      # Precision matrix
      x[[g]][[sigma]] <- as(solve_symmetric(x[[g]][[kappa]]),"Matrix")
      
      if (all) {
        x[[g]][[omega]] <- as(qgraph::wi2net(as.matrix(x[[g]][[kappa]])),"Matrix")
      }
    }
    
  }
  return(x)
}