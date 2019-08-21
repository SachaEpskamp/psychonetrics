matrixsetup_rho <- function(
  rho, # rho argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "rho",
  beta = array(0, c(nNode, nNode,nGroup))
){
  # Check if rho is character:
  ischar <- is.character(rho)
  
  # Fix rho:
  rho <- fixAdj(rho,nGroup,nNode,equal)
  
  # For each group, form starting values:
  rhoStart <- rho
  for (g in 1:nGroup){
    # Current estimate:
    corest <- cov2cor(as.matrix(spectralshift(expcov[[g]])))
    
    # Covs with sample covs:
    rhoStart[,,g] <-  1*(rhoStart[,,g]!=0) * corest  
    
    # If rho was a character, remove offdiagonal for endogenous variables:
    if (ischar && nNode > 1){
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(rho[,,g]) %in% endo | col(rho[,,g]) %in% endo) & (row(rho[,,g] ) != col(rho[,,g] ))
      rho[,,g][inds] <- rhoStart[,,g][inds] <-  0
    }
  }
  
  # Form the model matrix part:
  list(rho,
       mat =  name,
       op =  "~~",
       symmetrical= TRUE, 
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       diag0=TRUE,
       lower = -1,
       upper = 1,
       start = rhoStart,
       sampletable=sampletable
  )
}