matrixsetup_isingbeta <- function(
  beta, # sigma argument
  nGroup, # Number of groups
  equal = FALSE,
  sampletable,
  name = "beta"
){
  # Fix mu
  beta <- fixMu(beta,nGroup,1,equal)

  # For each group, form starting values:
  betaStart <- beta
  for (g in 1:nGroup){
    betaStart[,g][] <- 1
  }

  # Form the model matrix part:
  list(beta,
       mat =  name,
       op =  "temperature",
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = "beta",
       colnames = "",
       start = betaStart)
}