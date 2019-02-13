# Fix adjacency matrix structures:
fixAdj <- function(adjacency, nGroup, nNode, equal = FALSE){
  # Check adjacency:
  # Make adjacency if one of the defaults is used:
  if (is.character(adjacency)){
    
    # Empty network:
    if (adjacency == "empty"){
      
      # Equal for all groups:
      if (equal){
        adjacency <- array(diag(1+seq_len(nNode)), c(nNode, nNode, nGroup))
      } else {
        # Different for all groups:
        adjacency <- array(diag(nNode), c(nNode, nNode, nGroup))
      }
    } else {
      # Full network:
      if (equal){
        adjacency <- array(0, c(nNode, nNode, nGroup))
        for (i in 1:nGroup){
          adjacency[,,i][lower.tri(adjacency[,,i],diag=TRUE)] <- 1 + seq_len(sum(lower.tri(adjacency[,,i],diag=TRUE)))
          # adjacency[,,i][upper.tri(adjacency[,,i])] <- t(adjacency)[upper.tri(adjacency[,,i])] # Actually let's just ignore the upper.tri part!
        }
      } else {
        # Different for all groups:
        adjacency <- array(1, c(nNode, nNode, nGroup))
      }     
    }
  }
  
  
  # Check if the adjacency is a matrix:
  if (length(dim(adjacency)) == 2){
    adjacency <- array(adjacency, c(nNode, nNode, nGroup))
  }
  
  # Check dimensions:
  if (dim(adjacency)[1]!=nNode){
    stop("Number of rows in 'adjacency' does not equal the number of variables")
  }
  if (dim(adjacency)[2]!=nNode){
    stop("Number of columns in 'adjacency' does not equal the number of variables")
  }
  if (dim(adjacency)[3]!=nGroup){
    stop("Number of layers in 'adjacency' does not equal the number of groups")
  }
  
  # Clear al upper tris:
  for (i in 1:nGroup){
    adjacency[,,i][upper.tri(adjacency[,,i])] <- 0
  }
  
  # Should there be equality constrains?
  if (equal && any(adjacency==1)){
    curMax <- max(adjacency)
    adjacency[,,1][adjacency[,,1]==1&lower.tri(adjacency[,,1],diag=TRUE)] <- curMax + seq_len(sum(adjacency[,,1]==1&lower.tri(adjacency[,,1],diag=TRUE)))
    if (nGroup > 1){
      for (i in 2:nGroup){
        adjacency[,,i][adjacency[,,i]==1] <- adjacency[,,1][adjacency[,,i]==1]
      }      
    }
  }
  adjacency
}