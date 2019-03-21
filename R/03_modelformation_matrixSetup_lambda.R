matrixsetup_lambda <- function(
  lambda, # lambda argument
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  observednames, latentnames,
  equal = FALSE,
  sampletable,
  name = "lambda"
){
  # Fix lambda:
  lambda <- fixMatrix(lambda,nGroup,equal)
  
  # For each group, form starting values:
  lambdaStart <- lambda
  for (g in 1:nGroup){
    # Current cov estimate:
    curcov <- as.matrix(spectralshift(expcov[[g]]))
    
    # For all laodings:
    for (f in seq_len(ncol(lambdaStart))){
      # First principal component of sub cov:
      if (any(lambda[,f,g]!=0)){
        ev1 <- eigen(curcov[lambda[,f,g]!=0,lambda[,f,g]!=0])$vectors[,1]
        lambdaStart[lambdaStart[,f,g]!=0,f,g] <- ev1 / ev1[1]        
      } 
    }
  }
  
  # Form the model matrix part:
  list(lambda,
       mat = name,
       op =  "~=",
       rownames = observednames,
       colnames = latentnames,
       sparse = TRUE,
       start = lambdaStart,
       sampletable=sampletable
  )
}