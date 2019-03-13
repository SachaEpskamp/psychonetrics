matrixsetup_kappa <- function(
  kappa, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "kappa"
){
  # Fix lower tri:
  kappa <- fixAdj(kappa,nGroup,nNode,equal,diag0=FALSE)
  
  
  # For each group, form starting values:
  kappaStart <- kappa
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(spectralshift(expcov[[g]]))
    
    zeroes <- which(kappaStart[,,g]==0 & t(kappaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      wi <- corpcor::pseudoinverse(covest)
    } else {
      glas <- glasso(as.matrix(covest),
                     rho = 1e-10, zero = zeroes)
      wi <- glas$wi
    }
    
    # Network starting values:
    kappaStart[,,g] <- (kappaStart[,,g] !=0) * as.matrix(wi)
  }

  # Form the model matrix part:
  list(kappa,
       mat =  name,
       op =  "--",
       symmetrical= TRUE, 
       sampletable=sampletable,
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       diag0 = FALSE,
       start = kappaStart
  )
}