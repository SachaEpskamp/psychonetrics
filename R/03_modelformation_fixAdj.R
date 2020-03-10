# Fix kappa matrix structures:
fixAdj <- function(kappa, nGroup, nNode, equal = FALSE, diag0 = FALSE, diagonal = FALSE){
  # Check kappa:
  # Make kappa if one of the defaults is used:
  if (is.character(kappa)){
    
    # Empty network:
    if (kappa %in% c("diag","empty")){
      
      # Equal for all groups:
      if (equal){
        kappa <- array(diag(1+seq_len(nNode)), c(nNode, nNode, nGroup))
      } else {
        # Different for all groups:
        kappa <- array(diag(nNode), c(nNode, nNode, nGroup))
      }
    } else if (kappa == "zero"){
      
      kappa <- array(0, c(nNode, nNode, nGroup))
      
    } else {
      # Full network:
      if (equal){
        kappa <- array(0, c(nNode, nNode, nGroup))
        for (i in 1:nGroup){
          kappa[,,i][lower.tri(kappa[,,i],diag=TRUE)] <- 1 + seq_len(sum(lower.tri(kappa[,,i],diag=TRUE)))
          # kappa[,,i][upper.tri(kappa[,,i])] <- t(kappa)[upper.tri(kappa[,,i])] # Actually let's just ignore the upper.tri part!
        }
      } else {
        # Different for all groups:
        kappa <- array(1, c(nNode, nNode, nGroup))
      }     
    }
  }
  
  
  # Check if the kappa is a matrix:
  if (length(dim(kappa)) == 2){
    kappa <- array(kappa, c(nNode, nNode, nGroup))
  }
  
  # Check dimensions:
  if (dim(kappa)[1]!=nNode){
    stop("Number of rows in 'kappa' does not equal the number of variables")
  }
  if (dim(kappa)[2]!=nNode){
    stop("Number of columns in 'kappa' does not equal the number of variables")
  }
  if (dim(kappa)[3]!=nGroup){
    stop("Number of layers in 'kappa' does not equal the number of groups")
  }
  
  # Clear al upper tris and diag if needed:
  for (i in 1:nGroup){
    # Remove diagonal:
    if (diag0){
      kappa[,,i][upper.tri(kappa[,,i],diag=TRUE)] <- 0  
    } else if (diagonal){
      # Keep only diagonal:
      kappa[,,i][diag(nrow(kappa[,,i]))!=1] <- 0  
    } else {
      kappa[,,i][upper.tri(kappa[,,i])] <- 0
    }
  }
  
  # Should there be equality constrains?
  if (equal && any(kappa==1)){
    curMax <- max(kappa)
    kappa[,,1][kappa[,,1]==1&lower.tri(kappa[,,1],diag=TRUE)] <- curMax + seq_len(sum(kappa[,,1]==1&lower.tri(kappa[,,1],diag=TRUE)))
    if (nGroup > 1){
      for (i in 2:nGroup){
        kappa[,,i][kappa[,,i]==1] <- kappa[,,1][kappa[,,i]==1]
      }      
    }
  }
  kappa
}