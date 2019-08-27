# General matrix fixer
fixTau <- function(Matrix,sampleThresholds, equal = FALSE){

  nGroup <- length(sampleThresholds)
  if (length(unique(sapply(sampleThresholds,length))) > 1){
    stop("Different number of ordinal variables is not allowed.")
  }
  ncols <- length(sampleThresholds[[1]])
  
  # Thresholds per node:
  threshPerNode <- do.call(cbind,lapply(sampleThresholds,sapply,length))
  
  # Max number:
  nrows <- max(threshPerNode)

  # Check Matrix:
  # Make Matrix if one of the defaults is used:
  if (missing(Matrix)){

      # Full network:
      if (equal){
        Matrix <- array(0, c(nrows, ncols, nGroup))
        for (i in 1:nGroup){
          Matrix[,,i] <- 1 + seq_len(nrows*ncols)
        }
      } else {
        # Different for all groups:
        Matrix <- array(1, c(nrows, ncols, nGroup))
      }     
    
  } else {
    nrows <- dim(Matrix)[1]
    ncols <- dim(Matrix)[2]
  }
  
  # Check if the Matrix is a matrix:
  if (length(dim(Matrix)) == 2){
    Matrix <- array(Matrix, c(nrows, ncols, nGroup))
  }
  
  # Make the appropriate elements 0:
  for (g in seq_len(nGroup)){
    for (i in seq_len(ncols)){
      if (threshPerNode[i,g] < nrows){
        Matrix[(threshPerNode[i,g]+1):nrows,i,g] <- NA
      }

    }
  }

  # Check dimensions:
  if (dim(Matrix)[1]!=nrows){
    stop("Number of rows in 'Matrix' does not equal the number of variables")
  }
  if (dim(Matrix)[2]!=ncols){
    stop("Number of columns in 'Matrix' does not equal the number of variables")
  }
  if (dim(Matrix)[3]!=nGroup){
    stop("Number of layers in 'Matrix' does not equal the number of groups")
  }
  
  # # Clear al upper tris and diag if needed:
  # for (i in 1:nGroup){
  #   # Remove diagonal:
  #   if (diag0){
  #     Matrix[,,i][upper.tri(Matrix[,,i],diag=TRUE)] <- 0  
  #   } else if (diagonal){
  #     # Keep only diagonal:
  #     Matrix[,,i][diag(nrow(Matrix[,,i]))!=1] <- 0  
  #   } else {
  #     Matrix[,,i][upper.tri(Matrix[,,i])] <- 0
  #   }
  # }
  
  # Should there be equality constrains?
  if (equal && any(Matrix==1)){
    curMax <- max(Matrix)
    Matrix[,,1][Matrix[,,1]==1] <- curMax + seq_len(sum(Matrix[,,1]==1))
    if (nGroup > 1){
      for (i in 2:nGroup){
        Matrix[,,i][Matrix[,,i]==1] <- Matrix[,,1][Matrix[,,i]==1]
      }      
    }
  }
  Matrix
}