# Add standard errors and p-values!
# Add the modification indices:
addSEs <-  function(x){
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

  # Check if gradient and hessian are present:
  gradient <- !is.null(x@fitfunctions$gradient)
  hessian <- !is.null(x@fitfunctions$hessian)
  
  # Sample size:
  n <- sum(x@sample@groups$nobs)
  
  # Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but cosntrained per group #
  
  # # Compute a gradient:
  # if (gradient){
  #   g <- x@fitfunctions$gradient(parVector(modCopy), modCopy)
  # } else {
  #  g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  
  # Compute a hessian:
  if (hessian){
    H <- x@fitfunctions$hessian(parVector(x), x)
  } else if (gradient){
    H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(x), model=x) 
  } else {
    H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(x), model=x) 
  }
  
  # Obtain SEs
  # SEs <-  sqrt(abs(diag(solve(-n/2*H))))
  SEs <-  sqrt(abs(diag(solve(H))))
  # 
  # 
  # x <- sqrt(abs(diag(solve(H))))
  # sqrt(2/n) * x -
  # SEs
  # Add standard errors:
  for (i in unique(x@parameters$par[x@parameters$par!=0])){
    # Compute effective N:
    groups <- unique(x@parameters$group_id[x@parameters$par == i])
    Neff <- sum(x@sample@groups$nobs[x@sample@groups$id %in% groups])
    x@parameters$se[x@parameters$par == i] <- sqrt(2/Neff) * SEs[i]
  }
  
  # Add p-value:
  x@parameters$p <- 2*pnorm(abs(x@parameters$est),mean = 0,sd=x@parameters$se,lower.tail=FALSE)
  
  # FIXME: tempory round:
  # x@parameters$p <- round(x@parameters$p,5)
  
  # Return:
  return(x)
}