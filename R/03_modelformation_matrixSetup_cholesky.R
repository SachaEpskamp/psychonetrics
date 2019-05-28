matrixsetup_lowertri <- function(
  lowertri, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "lowertri",
  beta = array(0, c(nNode, nNode,nGroup))
){
  # Check if sigma is character:
  ischar <- is.character(lowertri)
  
  # Fix lower tri:
  lowertri <- fixAdj(lowertri,nGroup,nNode,equal)
  
  
  # For each group, form starting values:
  lowertriStart <- lowertri
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(expcov[[g]])
    
    # Start values:
    if (!any(is.na(covest))){
      tryres <- try(
        Lest <- t(as.matrix(chol(covest))), silent = TRUE
      )
      if (is(tryres, "try-error")){
        Lest <- diag(nrow(covest))
      }
      # Chol with sample cholesky:
      lowertriStart[,,g] <-  1*(lowertriStart[,,g]!=0) * Lest
    } else {
      lowertriStart[,,g] <-  1*(lowertriStart[,,g]!=0) * 0.05
    }
    
    if (ischar && nNode > 1){
      
      
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(lowertri[,,g]) %in% endo | col(lowertri[,,g]) %in% endo) & (row(lowertri[,,g] ) != col(lowertri[,,g] ))
      lowertri[,,g][inds] <- lowertriStart[,,g][inds] <-  0
    }
  }
  
  # Form the model matrix part:
  list(lowertri,
       mat =  name,
       op =  "~chol~",
       lowertri= TRUE, 
       sampletable=sampletable,
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       start = lowertriStart
  )
}