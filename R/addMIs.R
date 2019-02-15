# Add the modification indices:
addMIs <- addModificationIndices <- function(x){
  # If no constrained parameters, nothing to do!
  if (!any(x@parameters$par == 0)){
    return(x)
  }
  
  # Check if gradient and hessian are present:
  gradient <- !is.null(x@fitfunctions$gradient)
  hessian <- !is.null(x@fitfunctions$hessian)
  
  # Sample size:
  n <- sum(x@sample@groups$nobs)
  
  # Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but cosntrained per group #
  # Fully free:
  # Copy the model:
  modCopy <- x
  # Add free parameter numbers:
  modCopy@parameters$par[modCopy@parameters$par==0] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$par==0))
  
  # Remake the model matrix:
  modCopy@fitfunctions$extramatrices$M <- Mmatrix(modCopy@parameters)
  # Compute a gradient:

  if (gradient){
    g <- x@fitfunctions$gradient(parVector(modCopy), modCopy)
  } else {
   g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  }
  
  # Compute a hessian:
  if (hessian){
    H <- x@fitfunctions$hessian(parVector(modCopy), modCopy)
  } else if (gradient){
    H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(modCopy), model=modCopy) 
  } else {
    H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  }
  
  # For every new parameter:
  x@parameters$mi[] <- 0
  curMax <- max(x@parameters$par)
  for (i in which(x@parameters$par==0)){
    ind <- modCopy@parameters$par[i]
    curInds <- seq_len(curMax)
    x@parameters$mi[i] <- n * (0.5 * g[i]^2)/(H[ind,ind] - H[ind,curInds,drop=FALSE] %*% solve(H[curInds,curInds]) %*% H[curInds,ind,drop=FALSE])
    x@parameters$pmi[i] <- pchisq(x@parameters$mi[i],df = 1,lower.tail = FALSE)
  }
  
  # TODO CONSTRAINED MIs
  
  return(x)
}