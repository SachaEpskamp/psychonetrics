matrixsetup_sigma <- function(
  sigma, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "sigma",
  beta = array(0, c(nNode, nNode,nGroup))
){
  # FIXME: lowerbound for diagonal always 0:
  diaglower <- 0 # Lowerbound of the diagonal values
  
  
  # FIXME: correlations between endogenous latents are now removed even when sigma = "full"
  
  # Check if sigma is character:
  ischar <- is.character(sigma)
  
  # Fix Sigma:
  sigma <- fixAdj(sigma,nGroup,nNode,equal)
  
  # For each group, form starting values:
  lower <- sigmaStart <- sigma
  lower[] <- -Inf
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(spectralshift(expcov[[g]]))
    
    # Covs with sample covs:
    sigmaStart[,,g] <-  1*(sigmaStart[,,g]!=0) * covest  
    
    # Set lower bound:
    if (dim(lower)[1]>1){
      diag(lower[,,g]) <- diaglower  
    } else {
      lower[,,g] <- diaglower
    }
    
    
    # If Sigma was a character, remove off-diagonal for endogenous variables:
    if (ischar && nNode > 1){
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(sigma[,,g]) %in% endo | col(sigma[,,g]) %in% endo) & (row(sigma[,,g] ) != col(sigma[,,g] ))
      sigma[,,g][inds] <- sigmaStart[,,g][inds] <-  0
    }
  }
  
  # Form the model matrix part:
  list(sigma,
       mat =  name,
       op =  "~~",
       symmetrical= TRUE, 
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       start = sigmaStart,
       sampletable=sampletable,
       lower=lower
  )
}