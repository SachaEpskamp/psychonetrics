# Fix kappa matrix structures:
fixMu <- function(mu, nGroup, nNode, equal = FALSE){
  
  # Missing:
  if (missing(mu)){
    if (equal){
      mu <- matrix(1+seq_len(nNode),nNode,nGroup)
    } else {
      mu <- matrix(1,nNode,nGroup)      
    }
  }
  
  # If vector, repeat:
  if (!is.matrix(mu)){
    if (equal){
      mu <- matrix(mu,nNode,nGroup)
      curMax <- max(mu)
      mu[,1][mu[,1]==1] <- curMax + seq_len(sum(mu[,1]==1))
      if (nGroup > 1){
        for (i in 2:nGroup){
          mu[,i][mu[,i]==1] <- mu[,1][mu[,i]==1]  
        }
      }
    } else {
      mu <- matrix(mu,nNode,nGroup)
    }
  }
  
  # Check dimensions:
  if (!is.matrix(mu)){
    stop("'mu' must be a matrix")
  }
  if (nrow(mu) != nNode){
    stop("'mu' must have a row for each variable")
  }
  if (ncol(mu) != nGroup){
    stop("'mu' must have a column for each group")
  }
  
  # Return:
  mu
}