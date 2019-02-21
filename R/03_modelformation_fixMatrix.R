# General matrix fixer
fixMatrix <- function(Matrix, nGroup, nrows, ncols, equal = FALSE, diag0 = FALSE, diagonal = FALSE){
  if (missing(nrows) & is.character(Matrix)){
    stop("'nrows' must be assigned if Matrix is a character")
  }
  if (missing(ncols) & is.character(Matrix)){
    stop("'ncols' must be assigned if Matrix is a character")
  }

  # Check Matrix:
  # Make Matrix if one of the defaults is used:
  if (is.character(Matrix)){
    
    # Empty network:
    if (Matrix == "empty"){
        # Different for all groups:
        Matrix <- array(0, c(nrows, ncols, nGroup))
    } else {
      # Full network:
      if (equal){
        Matrix <- array(0, c(nrows, ncols, nGroup))
        for (i in 1:nGroup){
          Matrix[,,i] <- 1 + seq_len(nrows*ncols)
          # Matrix[,,i][upper.tri(Matrix[,,i])] <- t(Matrix)[upper.tri(Matrix[,,i])] # Actually let's just ignore the upper.tri part!
        }
      } else {
        # Different for all groups:
        Matrix <- array(1, c(nrows, ncols, nGroup))
      }     
    }
  } else {
    nrows <- dim(Matrix)[1]
    ncols <- dim(Matrix)[2]
  }
  
  
  # Check if the Matrix is a matrix:
  if (length(dim(Matrix)) == 2){
    Matrix <- array(Matrix, c(nrows, ncols, nGroup))
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