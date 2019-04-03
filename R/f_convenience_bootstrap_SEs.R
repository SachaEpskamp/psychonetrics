# Inner function:
bootstrap_SEs <- function(x,nCores = 1, reps = 1000, verbose = TRUE){
 pars <- x %>% replicator(bootstrap %>% runmodel, reps = reps, nCores = nCores, verbose = verbose, results = "parameters")
 x@parameters$se_boot <- ifelse(x@parameters$fixed,NA,sapply(pars,sd,na.rm=TRUE))
 x@parameters$p_boot <- ifelse(x@parameters$fixed,NA,2*pnorm(abs(x@parameters$est),mean = 0,sd=x@parameters$se_boot,lower.tail=FALSE))
 x
}