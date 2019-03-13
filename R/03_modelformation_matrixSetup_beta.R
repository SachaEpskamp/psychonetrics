matrixsetup_beta <- function(
  beta, # beta argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "beta"
){
  # Fix beta:
  beta <- fixMatrix(beta,nGroup = nGroup,nrows = nNode,ncols = nNode,equal = equal,diag0=TRUE)
  
  # For each group, form starting values:
  betaStart <- 0.1*(beta!=0)
  
  # Form the model matrix part:
  list(beta,
       mat =  name,
       op =  "<-",
       rownames = labels,
       colnames = labels,
       sparse = TRUE,
       start = betaStart,
       sampletable=sampletable
  )
}