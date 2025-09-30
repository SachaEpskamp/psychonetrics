matrixsetup_isingbeta <- function(
  beta, # sigma argument
  nGroup, # Number of groups
  equal = FALSE,
  sampletable,
  name,
  log = FALSE
){
  if (missing(name)){
     if (!log){
       name <- "beta"
     } else {
       name <- "log_beta"
     }
    
    
  }
  # Fix mu
  beta <- fixMu(beta,nGroup,1,equal)

  # For each group, form starting values:
  betaStart <- beta
  for (g in 1:nGroup){
    if (log){
      betaStart[,g][] <- 0  
    } else {
      betaStart[,g][] <- 1
    }
    
  }

  # Form the model matrix part:
  list(beta,
       mat =  name,
       op =  paste0(ifelse(log,"log ",""),"inverse temp"),
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = ifelse(log,"log_beta","beta"),
       colnames = "",
       start = betaStart)
}