matrixsetup_beta <- function(
  beta, # beta argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "beta",
  start
){
  # Fix beta:
  beta <- fixMatrix(beta,nGroup = nGroup,nrows = nNode,ncols = nNode,equal = equal,diag0=TRUE)
  
  # For each group, form starting values:
  if (missing(start)){
    betaStart <- 0.1*(beta!=0)  
  } else {
    betaStart <- beta
    for (g in seq_along(start)){
      betaStart[,,g] <- (beta[,,g]!=0) * start[[g]]
    }
  }
  
  
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