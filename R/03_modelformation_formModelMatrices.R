# General function that forms the model matrices
formModelMatrices <- function(x){
  # x is the model:
  pars <- x@parameters
  mats <- x@matrices
  
  # All groups:
  allGroups <- unique(pars$group_id)
  nGroup <- length(allGroups)

  # form matrices:
  Matrices <- lapply(seq_len(nGroup),function(g){
    groupMod <- lapply(seq_len(nrow(mats)),function(i){
      mat <- matrix(0, mats$nrow[i], mats$ncol[i])
      for (id in which(pars$matrix==mats$name[i]&pars$group_id==g)){
        mat[pars$row[id],pars$col[id]] <- pars$est[id]
      }
      if (mats$symmetrical[i]){
        mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
      } 
      # What kind of matrix?
      if (mats$diagonal[i]){
        mat <- as(mat, "diagonalMatrix")
        
      } else if (mats$lowertri[i]){
        mat <- as(mat, "Matrix")
      } else if (mats$sparse[i]){
        if (mats$symmetrical[i]){
          mat <- as(mat, "dsCMatrix")
        } else {
          mat <- as(mat, "sparseMatrix")
        }
      } else if (mats$posdef[i]){
       mat <- as(mat, "dpoMatrix")
      } else if (mats$symmetrical[i]){
      mat <- as(mat, "dsyMatrix")
      } else {
        mat <- as(mat, "dgeMatrix")
      }
      
      mat
    })
    names(groupMod) <-  mats$name
    groupMod
  })
  names(Matrices) <- allGroups
  
  # Some additional matrices need to be formed:
  if (x@model == "varcov"){
    # Variance covariance matrices:
    for (g in seq_len(nGroup)){
      if (x@types$y == "chol"){
        Matrices[[g]]$sigma <- as(Matrices[[g]]$lowertri %*% t(Matrices[[g]]$lowertri), "dsyMatrix")
        # Matrices[[g]]$rho <- cov2cor(Matrices[[g]]$sigma)
        Matrices[[g]]$kappa <- as(corpcor::pseudoinverse(Matrices[[g]]$sigma), "sparseMatrix")
        Matrices[[g]]$omega <- as(qgraph::wi2net(as.matrix(Matrices[[g]]$kappa)),"sparseMatrix")
      } else if (x@types$y == "cov"){
        # Matrices[[g]]$rho <- cov2cor(Matrices[[g]]$sigma)
        Matrices[[g]]$kappa <- as(corpcor::pseudoinverse(Matrices[[g]]$sigma), "sparseMatrix")
        Matrices[[g]]$omega <- as(qgraph::wi2net(as.matrix(Matrices[[g]]$kappa)),"sparseMatrix")
      } else if (x@types$y == "ggm"){
        Matrices[[g]]$sigma <- Matrices[[g]]$delta %*% corpcor::pseudoinverse(spectralshift(Diagonal(ncol(Matrices[[g]]$omega)) - Matrices[[g]]$omega)) %*% Matrices[[g]]$delta
        Matrices[[g]]$kappa <- corpcor::pseudoinverse(Matrices[[g]]$sigma)
      } else if (x@types$y == "prec"){
        Matrices[[g]]$sigma <- corpcor::pseudoinverse(spectralshift(Matrices[[g]]$kappa))
        Matrices[[g]]$omega <- as(qgraph::wi2net(as.matrix(Matrices[[g]]$kappa)),"sparseMatrix")
      }
    }
  }

  return(Matrices)
}