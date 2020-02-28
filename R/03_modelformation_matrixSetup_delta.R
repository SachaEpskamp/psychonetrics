matrixsetup_delta <- function(
  delta, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "delta",
  onlyStartSign = FALSE,
  omegaStart
){
  # Fix lower tri:
  delta <- fixAdj(delta,nGroup,nNode,equal,diagonal=TRUE)
  
  
  # For each group, form starting values:
  deltaStart <- delta
  
  # OmegaStart:
  if (missing(omegaStart)){
    omegaStart <- array(1,c(nNode,nNode,nGroup))
  }
  
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(spectralshift(expcov[[g]]))
    
    # FIXME: This is already done in omega computation, but I am lazy....
    zeroes <- which(omegaStart[,,g]==0 & t(omegaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      wi <- solve_symmetric(spectralshift(covest))
    } else {
      glas <- glasso(as.matrix(covest),
                     rho = 1e-10, zero = zeroes)
      wi <- glas$wi
    }
    
    # Network starting values:
    deltaStart[,,g] <- diag(1/sqrt(diag(wi)))
    
    if (onlyStartSign){
      deltaStart[,,g] <- ifelse(deltaStart[,,g]!=0,sign(deltaStart[,,g]), 0)
    }
  }
  
  
  # Form the model matrix part:

  list(delta,
       mat =  name,
       op =  "~/~",
       symmetrical= TRUE, 
       sampletable=sampletable,
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       diagonal = TRUE,
       lower = 0,
       start = deltaStart
  )
}