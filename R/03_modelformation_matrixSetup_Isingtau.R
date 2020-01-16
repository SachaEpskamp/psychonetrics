matrixsetup_isingtau <- function(
  tau, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "tau"
){
  # Fix mu
  tau <- fixMu(tau,nGroup,nNode,equal)

  # For each group, form starting values:
  tauStart <- tau
  for (g in 1:nGroup){
    tauStart[,g][] <- 0
  }

  # Form the model matrix part:
  list(tau,
       mat =  name,
       op =  "~1",
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = labels,
       colnames = "1",
       start = tauStart)
}