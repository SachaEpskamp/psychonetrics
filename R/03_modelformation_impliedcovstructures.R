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
    
    rho = paste0("rho","_",name)
    SD = paste0("SD","_",name)
    
    IplusRho = paste0("IplusRho","_",name)
    SD_IplusRho = paste0("SD_IplusRho","_",name)
    
    
  } else {
    sigma = "sigma"
    omega = "omega"
    delta = "delta"
    kappa = "kappa"
    lowertri = "lowertri"
    IminOinv = "IminOinv"
    delta_IminOinv = "delta_IminOinv"
  
    rho = "rho"
    SD = "SD"
    IplusRho = "IplusRho"
    SD_IplusRho = "SD_IplusRho"
  }
  
  for (g in seq_len(nGroup)){
    # Form the models:
    # Contemporaneous:
    if (type == "cov"){
      
      # Only need to do things if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- sparseordense(solve_symmetric(x[[g]][[sigma]]))
          x[[g]][[omega]]  <- sparseordense(qgraph::wi2net(as.matrix(x[[g]][[kappa]])))
          x[[g]][[rho]] <- sparseordense(cov2cor(as.matrix(x[[g]][[sigma]])))
          x[[g]][[SD]] <- Diagonal(x=diag(sqrt(diag(as.matrix(x[[g]][[sigma]])))))
        }
      }
    } else if(type == "chol"){
      # form cov matrix:
      x[[g]][[sigma]] <- sparseordense(x[[g]][[lowertri]] %*% t(x[[g]][[lowertri]]))
      
      # Return precision and network if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- sparseordense(solve_symmetric(x[[g]][[sigma]]))
          x[[g]][[omega]]  <- sparseordense(qgraph::wi2net(as.matrix(x[[g]][[kappa]])))
          x[[g]][[rho]] <- sparseordense(cov2cor(as.matrix(x[[g]][[sigma]])))
          x[[g]][[SD]] <- Diagonal(x=diag(sqrt(diag(as.matrix(x[[g]][[sigma]])))))
        }

      }
    } else if (type == "ggm"){
      # First check if the delta Matrix is present (it is ignored when corinput = TRUE only, so don't need to know that that argument was used):
      if (is.null(x[[g]][[delta]])){
        # FIXME: non positive definite matrices are even worse here... So I am trying to solve this with a spectral shift for now:
        IminO_dummy <-  sparseordense(solve_symmetric(spectralshift(Diagonal(ncol(x[[g]][[omega]])) - x[[g]][[omega]])))
        x[[g]][[delta]] <- Diagonal(x = diag(IminO_dummy)^(-0.5))
      } else {
        IminO_dummy <-  sparseordense(solve_symmetric(Diagonal(ncol(x[[g]][[omega]])) - x[[g]][[omega]]))
        
      }
      
      x[[g]][[sigma]] <- sparseordense(x[[g]][[delta]] %*%IminO_dummy  %*% x[[g]][[delta]])
      
      # Stuff needed if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- sparseordense(solve_symmetric(x[[g]][[sigma]]))
          x[[g]][[rho]] <- sparseordense(cov2cor(as.matrix(x[[g]][[sigma]])))
          x[[g]][[SD]] <- Diagonal(x=diag(sqrt(diag(as.matrix(x[[g]][[sigma]])))))
        }
        
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]][[IminOinv]] <- IminO_dummy
        x[[g]][[delta_IminOinv]] <- sparseordense(x[[g]][[delta]] %*% x[[g]][[IminOinv]])
      }
    } else if (type == "prec"){
      # Precision matrix
      x[[g]][[sigma]] <- sparseordense(solve_symmetric(x[[g]][[kappa]]))
      
      if (all) {
        x[[g]][[omega]] <- sparseordense(qgraph::wi2net(as.matrix(x[[g]][[kappa]])))
        x[[g]][[rho]] <- sparseordense(cov2cor(as.matrix(x[[g]][[sigma]])))
        x[[g]][[SD]] <- Diagonal(x=diag(sqrt(diag(as.matrix(x[[g]][[sigma]])))))
      }
    } else if (type == "cor"){

      # First check if the SD Matrix is present (it is ignored when corinput = TRUE only, so don't need to know that that argument was used):
      if (is.null(x[[g]][[SD]])){
        x[[g]][[SD]] <- Diagonal(ncol(x[[g]][[rho]]))
      }
      
      x[[g]][[sigma]] <- sparseordense(x[[g]][[SD]] %*% (Diagonal(ncol(x[[g]][[rho]])) + x[[g]][[rho]]) %*% x[[g]][[SD]])
      
      # Stuff needed if all = TRUE:
      if (all){
        if (!all(x[[g]][[sigma]] == 0)){
          x[[g]][[kappa]] <- sparseordense(solve_symmetric(x[[g]][[sigma]]))
          x[[g]][[omega]]  <- sparseordense(qgraph::wi2net(as.matrix(x[[g]][[kappa]])))
        }
        
      }
      
      # Extra matrix needed:
      if (!all){
        x[[g]][[IplusRho]] <- sparseordense(Diagonal(ncol(x[[g]][[rho]])) + x[[g]][[rho]])
        x[[g]][[SD_IplusRho]] <- sparseordense(x[[g]][[SD]] %*% x[[g]][[IplusRho]])
      }
      
      
      
    }
    
  }
  return(x)
}
