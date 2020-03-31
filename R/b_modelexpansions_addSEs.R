# Add standard errors and p-values!
# Add the modification indices:
addSEs <-  function(x,
                    verbose
                    # approxHessian = FALSE # Approximate the hessian even if it is already stored
                    ){
  
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  if (verbose){
    message("Adding standard errors...")
  }
  
  # If not computed, warn user!
  if (!x@computed){
    warning("Model was not computed, interpret standard errors and p-values with care!")
  }
  
  # If no free parameters, nothing to do!
  if (all(x@parameters$par == 0)){
    return(x)
  }
  

  # Clear old:
  x@parameters$se[] <- NA
  x@parameters$p[] <- NA

  # Sample size:
  n <- sum(x@sample@groups$nobs)
  
  # Compute a hessian if requested or needed:
  # if (!is.null(x@information)){
  #   Hinv <- solve_symmetric(x@information)
  # } else {
  #   Hinv <- solve_symmetric(psychonetrics_FisherInformation(x))
  # } 
  Hinv <- getVCOV(x)
  # if (!is.null(x@information)){
  #   Hinv <- solve_symmetric(x@information)
  # } else if (!is.null(x@optim$inverseHessian)){
  #   Hinv <- x@optim$inverseHessian
  # } else {
  #   # Check if gradient and hessian are present:
  #   gradient <- !is.null(x@fitfunctions$gradient)
  #   hessian <- !is.null(x@fitfunctions$hessian)
  #   
  #   if (hessian){
  #     H <- x@fitfunctions$hessian(parVector(x), x)
  #   } else if (gradient){
  #     H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(x), model=x) 
  #   } else {
  #     H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(x), model=x) 
  #   }
  #   
  #   # Invert:
  #   Hinv <- solve_symmetric(H)
  # } 
 
  # Obtain SEs
  # SEs <-  sqrt(abs(diag(solve_symmetric(-n/2*H))))
  # Hinv <- solve_symmetric(x@fitfunctions$information(x))
  SEs <-  sqrt(abs(diag(Hinv)))
  
  # 
  # 
  # x <- sqrt(abs(diag(solve_symmetric(H))))
  # sqrt(2/n) * x -
  # SEs
  N <- sum(x@sample@groups$nobs)
  # Add standard errors:
  for (i in unique(x@parameters$par[x@parameters$par!=0])){
    # Compute effective N:
    groups <- unique(x@parameters$group_id[x@parameters$par == i])
    Neff <- sum(x@sample@groups$nobs[x@sample@groups$id %in% groups])
    # x@parameters$se[x@parameters$par == i] <- sqrt(2/(Neff/N)) * SEs[i]
    x@parameters$se[x@parameters$par == i] <- SEs[i]
  }
  
  # Add p-value:
  x@parameters$p <- 2*pnorm(abs(x@parameters$est),mean = 0,sd=x@parameters$se,lower.tail=FALSE)
  
  # FIXME: tempory round:
  # x@parameters$p <- round(x@parameters$p,5)
  
  # Return:
  return(x)
}