matrixsetup_sigma <- function(
  sigma, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable
){
  # Fix Sigma:
  sigma <- fixAdj(sigma,nGroup,nNode,equal)
  
  # For each group, form starting values:
  sigmaStart <- sigma
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(spectralshift(expcov[[g]]))
    
    # Covs with sample covs:
    sigmaStart[,,g] <-  1*(sigmaStart[,,g]!=0) * covest  
  }
  
  # Form the model matrix part:
  list(sigma,
       mat =  "sigma",
       op =  "~~",
       symmetrical= TRUE, 
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       start = sigmaStart,
       sampletable=sampletable
  )
}