# General function that forms the model matrices
formModelMatrices <- function(x){
  # x is the model:
  pars <- x@parameters
  mats <- x@matrices
  
  # All groups:
  allGroups <- x@sample@groups$id
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
      # if (mats$diagonal[i]){
      #   mat <- as(mat, "diagonalMatrix")
      #   
      # } else if (mats$lowertri[i]){
      #   mat <- as(mat, "Matrix")
      # } else if (mats$sparse[i]){
      #   if (mats$symmetrical[i]){
      #     mat <- as(mat, "dsCMatrix")
      #   } else {
      #     mat <- as(mat, "sparseMatrix")
      #   }
      # } else if (mats$posdef[i]){
      #  mat <- as(mat, "dpoMatrix")
      # } else if (mats$symmetrical[i]){
      # mat <- as(mat, "dsyMatrix")
      # } else {
      #   mat <- as(mat, "dgeMatrix")
      # }
      mat <- as(mat, "Matrix")
      
      mat
    })
    names(groupMod) <-  mats$name
    groupMod
  })

  names(Matrices) <-  unique(x@sample@groups$label)
  

  return(Matrices)
}