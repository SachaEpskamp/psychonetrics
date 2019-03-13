# Implied model for precision. Requires appropriate model matrices:
implied_varcov <- function(model, all = FALSE){
  x <- formModelMatrices(model)
  
  for (g in seq_along(x)){
    if (model@types$y == "chol"){
      x[[g]]$sigma <- as(x[[g]]$lowertri %*% t(x[[g]]$lowertri), "Matrix")
      # Matrices[[g]]$rho <- cov2cor(Matrices[[g]]$sigma)
      x[[g]]$kappa <- as(corpcor::pseudoinverse(x[[g]]$sigma), "Matrix")
      if (all) x[[g]]$omega <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa)),"sparseMatrix")
    } else if (model@types$y == "cov"){
      # Matrices[[g]]$rho <- cov2cor(Matrices[[g]]$sigma)
      x[[g]]$kappa <- as(corpcor::pseudoinverse(x[[g]]$sigma), "Matrix")
      if (all) x[[g]]$omega <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa)),"sparseMatrix")
    } else if (model@types$y == "ggm"){
      x[[g]]$sigma <- x[[g]]$delta %*% corpcor::pseudoinverse(spectralshift(Diagonal(ncol(x[[g]]$omega)) - x[[g]]$omega)) %*% x[[g]]$delta
      x[[g]]$kappa <- corpcor::pseudoinverse(x[[g]]$sigma)
      
      # One usefull matrix in computation, but ignored in final:
      if (!all) x[[g]]$IminOinv <- corpcor::pseudoinverse(Diagonal(nrow(x[[g]]$omega)) - x[[g]]$omega)
      
    } else if (model@types$y == "prec"){
      x[[g]]$sigma <- corpcor::pseudoinverse(spectralshift(x[[g]]$kappa))
      if (all) x[[g]]$omega <- as(qgraph::wi2net(as.matrix(x[[g]]$kappa)),"sparseMatrix")
    }
  }
  
  x
  
  #   # For each group:
  # nGroup <- length(x)
  # Result <- lapply(seq_len(nGroup), function(g){
  #   
  #   # Implied precision:
  #   kappa <- corpcor::pseudoinverse(x[[g]]$sigma)
  # 
  #   # Implied means
  #   # mu <- x[[g]]$mu
  #   
  #   return(list(
  #     kappa = kappa
  #     # mu = mu
  #   )
  #   )
  # })
  # 
  # Result
}
