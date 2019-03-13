matrixsetup_mu <- function(
  mu, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expmeans, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "mu"
){
  # Fix mu
  mu <- fixMu(mu,nGroup,nNode,equal)

  # For each group, form starting values:
  muStart <- mu
  for (g in 1:nGroup){
    # Current estimate:
    meanest <-  as.matrix(expmeans[[g]])

    # Covs with sample covs:
    muStart[,g] <- 1*(mu[,g]!=0) * meanest
  }

  # Form the model matrix part:
  list(mu,
       mat =  name,
       op =  "~1",
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = labels,
       colnames = "1",
       start = muStart)
}