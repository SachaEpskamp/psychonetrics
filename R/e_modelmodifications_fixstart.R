# Function to fix start values:
fixstart <- function(x, reduce = 0.5, maxdiff = 0.1, tol =  0.01, maxtry = 25){
  
  stopifnot(is(x,"psychonetrics"))
  
  # Counter:
  gen <- 0
  check <- checkJacobian(x,plot = FALSE)
  
  repeat{
    if (mean(abs(check$numeric - check$analytic)) < 0.0001){
      break
    } else {
      gen <- gen + 1
      if (gen > maxtry ) stop("'maxtry' iteration reached without fixing start values.")
    }
    
    # Gradient recovery mechanism:
    
    freepars <- x@parameters[match(seq_len(max(x@parameters$par)),x@parameters$par),]
    freepars$diff <- abs(check$analytic - check$numeric)/abs(check$numeric)
    
    # Adjust the starting values of the parameters with largerst analytic-numeric differences:
    adjust <- freepars$par[freepars$diff>maxdiff & (
      grepl("beta",freepars$matrix) | (
      freepars$row != freepars$col & (
      grepl("omega",freepars$matrix)|grepl("lowertri",freepars$matrix)|grepl("sigma",freepars$matrix)|grepl("kappa",freepars$matrix)|grepl("rho",freepars$matrix)
      ))
    )]
    x@parameters$est[x@parameters$par %in% adjust] <- reduce * x@parameters$est[x@parameters$par %in% adjust]
    check <- checkJacobian(x,plot = FALSE)

  }
  
  x@modelmatrices <- formModelMatrices(x)
  return(x)
}