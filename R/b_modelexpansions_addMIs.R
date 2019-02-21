# Full function:
addMIs <- addModificationIndices <- function(x){
 x %>% addMIs_inner(equal=FALSE) %>% addMIs_inner(equal=TRUE) 
}

# Add the modification indices:
addMIs_inner <- addModificationIndices_inner <- function(x, equal = FALSE){
  # If no constrained parameters, nothing to do!
  if (!any(x@parameters$par == 0)){
    return(x)
  }
  
  # Clear old MIs:
  if (!equal){
    x@parameters$mi[] <- 0
    x@parameters$pmi[] <- NA
  } else {
    x@parameters$mi_equal[] <- 0
    x@parameters$pmi_equal[] <- NA
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

  # Obtain the full set of parameters that are constrained across all groups:
  if (equal){
    sum <- modCopy@parameters %>% group_by_("matrix","row","col") %>% summarize_(anyConstrained = ~any(fixed)) %>% 
      filter_(~anyConstrained)
    # Add a unique number to each:
    sum$par2 <-  max(modCopy@parameters$par) + seq_len(nrow(sum))
    
    # Left join back:
    modCopy@parameters <- modCopy@parameters %>% left_join(sum,by=c("matrix","row","col")) %>% 
      mutate(par = ifelse(par==0,par2,par))
    
  } else {
    # Add free parameter numbers:
    modCopy@parameters$par[modCopy@parameters$par==0] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$par==0))
    
    # For each group, free all parameters from equality constraints:
    if (nrow(modCopy@sample@groups)>1){
      for (g in 2:nrow(modCopy@sample@groups)){
        modCopy@parameters$par[modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par)] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par)))
      }
    }
  }

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
  curMax <- max(x@parameters$par)

  
  # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & x@parameters$par == 0]))){
  # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & (x@parameters$par == 0 | duplicated(x@parameters$par)|rev(duplicated(rev(x@parameters$par))))]))){
  for (i in sort(unique(modCopy@parameters$par))){
    ind <- i
    curInds <- seq_len(curMax)
    curInds <- curInds[curInds != i]
    mi <- n * (0.5 * g[i]^2)/(H[ind,ind] - H[ind,curInds,drop=FALSE] %*% solve(H[curInds,curInds]) %*% H[curInds,ind,drop=FALSE])
    p <- pchisq(mi,df = 1,lower.tail = FALSE)      
    if (equal){
      x@parameters$mi_equal[modCopy@parameters$par == i] <- round(mi,10) # round(mi,3)
      x@parameters$pmi_equal[modCopy@parameters$par == i] <- round(p, 10)
    } else {
      x@parameters$mi[modCopy@parameters$par == i] <- round(mi,10) # round(mi, 3)
      x@parameters$pmi[modCopy@parameters$par == i] <- round(p,10)
    }

  }
  
  if (equal){
    x@parameters$mi_equal[is.na(x@parameters$mi_equal)] <- 0
    x@parameters$pmi_equal[is.na(x@parameters$pmi_equal)] <- 1
  } else {
    x@parameters$mi[is.na(x@parameters$mi)] <- 0
    x@parameters$pmi[is.na(x@parameters$pmi)] <- 1
  }

  return(x)
}