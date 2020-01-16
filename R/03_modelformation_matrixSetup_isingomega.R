matrixsetup_isingomega <- function(
  omega, # omega argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "omega"
){
  # Check if omega is character:
  ischar <- is.character(omega)
  
  # Fix Omega:
  omega <- fixAdj(omega,nGroup,nNode,equal)
  
  # For each group, form starting values:
  omegaStart <- omega
  omegaStart[] <- 0
  
  # Form the model matrix part:
  list(omega,
       mat =  name,
       op =  "--",
       symmetrical= TRUE, 
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       posdef = FALSE,
       diag0=TRUE,
       start = omegaStart,
       sampletable=sampletable
  )
}