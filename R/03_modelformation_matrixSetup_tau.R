# Setup thresholds

matrixsetup_tau <- function(
  tau, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  sampleThresholds, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "tau"
){
  # Fix mu
  tau <- fixTau(tau,sampleThresholds,equal)

  # For each group, form starting values:
  tauStart <- tau
  for (g in 1:nGroup){
    # Current estimate:
    tauest <- matrix(0,nrow(tauStart),ncol(tauStart)) 
    for (i in 1:ncol(tau)){
      tauest[1:length(sampleThresholds[[g]][[i]]),i] <- sampleThresholds[[g]][[i]]
    }
    # Covs with sample covs:
    tauStart[,,g] <- 1*(tau[,,g]!=0) * tauest
  }

  # Form the model matrix part:
  list(tau,
       mat =  name,
       op =  "|",
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = paste0("t",seq_len(ncol(tau))),
       colnames = labels,
       start = tauStart,
       incomplete = TRUE)
}