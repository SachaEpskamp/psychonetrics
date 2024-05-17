matrixsetup_omega <- function(
    omega, # sigma argument
    nNode, # Number of nodes
    nGroup, # Number of groups
    expcov, # Expected covariance matrices (list).
    labels,
    equal = FALSE,
    sampletable,
    name = "omega",
    beta = array(0, c(nNode, nNode,nGroup)),
    onlyStartSign = FALSE,
    lassofix = TRUE 
){
  
  # Check if sigma is character:
  ischar <- is.character(omega)
  
  # Fix lower tri:
  omega <- fixAdj(omega,nGroup,nNode,equal,diag0=TRUE)
  
  
  # For each group, form starting values:
  omegaStart <- omega
  for (g in 1:nGroup){
    
    # If nNode == 1 just skip:
    if (nNode == 1){
      omegaStart[,,g] <- 0
      next
    }
    
    # Current estimate:
    covest <- as.matrix(expcov[[g]])
    
    
    if (!any(is.na(covest))){
      
      zeroes <- which(omegaStart[,,g]==0 & t(omegaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
      if (nrow(zeroes) == 0){
        wi <- solve_symmetric_cpp_matrixonly(covest)
        pcor <- qgraph::wi2net(as.matrix(wi))
        
        # FIXME: Quick check, if there is an outrageous starting value, use glasso with lasso instead:
        if (lassofix && (any(is.na(pcor)) || any(abs(pcor) > 0.8))){
          wi <- glasso(as.matrix(spectralshift(covest)), rho = 0.1)$wi
          pcor <-  qgraph::wi2net(as.matrix(wi))
        }
        
      } else {
        glas <- glasso(as.matrix(covest),
                       rho = 1e-10, zero = zeroes)
        wi <- glas$wi
        pcor <- qgraph::wi2net(as.matrix(wi))
      }
      
      # Network starting values:
      omegaStart[,,g] <- as.matrix(pcor)
      diag(omegaStart[,,g] ) <- 0  
      
      
      
      
    } else {
      
      # Network starting values:
      omegaStart[,,g] <- (omegaStart[,,g]!=0) * 0.001
      diag(omegaStart[,,g] ) <- 0
    }
    
    if (onlyStartSign){
      omegaStart[,,g] <- ifelse(omegaStart[,,g]!=0,0.001 * sign(omegaStart[,,g]), 0)
    }
    
    # If omega was a character, remove offdiagonal for endogenous variables:
    if (ischar && nNode > 1){
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(omega[,,g]) %in% endo | col(omega[,,g]) %in% endo) & (row(omega[,,g] ) != col(omega[,,g] ))
      omega[,,g][inds] <- omegaStart[,,g][inds] <-  0
    }
    
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