matrixsetup_SD <- function(
  SD, # SD argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "SD",
  beta = array(0, c(nNode, nNode,nGroup))
){
  # Check if SD is character:
  ischar <- is.character(SD)
  
  # Fix SD:
  SD <- fixAdj(SD,nGroup,nNode,equal)
  
  # For each group, form starting values:
  SDStart <- SD
  for (g in 1:nGroup){
    # Current estimate:
    SDest <- diag(sqrt(diag(as.matrix(spectralshift(expcov[[g]])))))
    
    # Covs with sample covs:
    SDStart[,,g] <-  1*(SDStart[,,g]!=0) * SDest  
    
    # If SD was a character, remove offdiagonal for endogenous variables:
    if (ischar && nNode > 1){
      # Which are endogenous?
      endo <- which(rowSums(beta[,,g])>0)
      
      # Remove these:
      inds <- (row(SD[,,g]) %in% endo | col(SD[,,g]) %in% endo) & (row(SD[,,g] ) != col(SD[,,g] ))
      SD[,,g][inds] <- SDStart[,,g][inds] <-  0
    }
  }
  
  # Form the model matrix part:
  list(SD,
       mat =  name,
       op =  "~~",
       symmetrical= TRUE, 
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = TRUE,
       diagonal = TRUE,
       lower = 0,
       start = SDStart,
       sampletable=sampletable
  )
}