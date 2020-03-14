matrixsetup_kappa <- function(
  kappa, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "kappa",
  beta = array(0, c(nNode, nNode,nGroup))
){
  # Check if kappa is character:
  ischar <- is.character(kappa)
  
  # Fix lower tri:
  kappa <- fixAdj(kappa,nGroup,nNode,equal,diag0=FALSE)
  
  
  # For each group, form starting values:
  kappaStart <- kappa
  for (g in 1:nGroup){
    # Current estimate:
    covest <- as.matrix(expcov[[g]])
    
    zeroes <- which(kappaStart[,,g]==0 & t(kappaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      wi <- solve_symmetric(covest)
      
      # FIXME: Quick check, if there is an outrageous starting value, use glasso with lasso instead:
      if (any(abs(qgraph::wi2net(wi)[lower.tri(wi,diag=FALSE)]) > 0.8)){
        wi <- glasso(as.matrix(spectralshift(covest)), rho = 0.1)$wi
      }
      
    } else {
      glas <- glasso(as.matrix(covest),
                     rho = 1e-10, zero = zeroes)
      wi <- glas$wi
    }
    
    # Network starting values:
    kappaStart[,,g] <- (kappaStart[,,g] !=0) * as.matrix(wi)
    
    # If kappa was a character, remove offdiagonal for endogenous variables:
    if (ischar && nNode > 1){
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(kappa[,,g]) %in% endo | col(kappa[,,g]) %in% endo) & (row(kappa[,,g] ) != col(kappa[,,g] ))
      kappa[,,g][inds] <- kappaStart[,,g][inds] <-  0
    }
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