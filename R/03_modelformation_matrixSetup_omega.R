matrixsetup_omega <- function(
  omega, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "omega"
){
  # Fix lower tri:
  omega <- fixAdj(omega,nGroup,nNode,equal,diag0=TRUE)
  
  
  # For each group, form starting values:
  omegaStart <- omega
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(spectralshift(expcov[[g]]))
    
    zeroes <- which(omegaStart[,,g]==0 & t(omegaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      wi <- corpcor::pseudoinverse(covest)
    } else {
      glas <- glasso(as.matrix(covest),
                     rho = 1e-10, zero = zeroes)
      wi <- glas$wi
    }
    
    # Network starting values:
    omegaStart[,,g] <- as.matrix(qgraph::wi2net(wi))
    diag(omegaStart[,,g] ) <- 0
  }
  
  # Form the model matrix part:
  list(omega,
       mat =  name,
       op =  "--",
       symmetrical= TRUE, 
       sampletable=sampletable,
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       diag0=TRUE,
       lower = -1,
       upper = 1,
       start = omegaStart
  )
}