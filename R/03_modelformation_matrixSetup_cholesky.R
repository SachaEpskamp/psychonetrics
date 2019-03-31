matrixsetup_lowertri <- function(
  lowertri, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "lowertri"
){
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