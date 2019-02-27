# Full function:
addMIs <- function(x, matrices = "all", type =  c("normal","free","equal"),verbose = TRUE){
  # full <- TRUE
  # if (full){
  if (verbose){
    message("Computing modification indices...")
  }
  
  if ("normal" %in% type){

    x <-  x %>% addMIs_inner_full(type = "normal")     
  }
  
  if (nrow(x@sample@groups) > 1){
    if ("free" %in% type){
      # if (verbose){
      #   message("Computing constrain-free modification indices...")
      # }
      x <- x %>% addMIs_inner_full(type = "free")         
    }
    if ("equal" %in% type){
      
      # if (verbose){
      #   message("Computing group-constrained indices...")
      # }
      x <- x %>% addMIs_inner_full(type = "equal")       
    }
    
  }
  return(x)

}


# Add the modification indices (FULL VERSION):
addMIs_inner_full <- function(x, type =  c("normal","free","equal")){
  type <- match.arg(type)
  
  # If no constrained parameters, nothing to do!
  if (!any(x@parameters$par == 0)){
    return(x)
  }
  
  # Clear old MIs:
  if (type == "normal"){
    x@parameters$mi[] <- 0
    x@parameters$pmi[] <- NA
    x@parameters$epc[] <- NA
  } else if (type == "free"){
    x@parameters$mi_free[] <- 0
    x@parameters$pmi_free[] <- NA
    x@parameters$epc_free[] <- NA
  } else {
    x@parameters$mi_equal[] <- 0
    x@parameters$pmi_equal[] <- NA
    x@parameters$epc_equal[] <- NA
  }
  
  # Sample size:
  n <- sum(x@sample@groups$nobs)
  
  # Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but constrained per group #
  # Fully free:
  # Copy the model:
  modCopy <- x
  
  # Obtain the full set of parameters that are constrained across all groups:
  if (type == "equal"){
    sum <- modCopy@parameters %>% group_by_("matrix","row","col") %>% summarize_(anyConstrained = ~any(fixed)) %>% 
      filter_(~anyConstrained)
    # Add a unique number to each:
    sum$par2 <-  max(modCopy@parameters$par) + seq_len(nrow(sum))
    
    # Left join back:
    modCopy@parameters <- modCopy@parameters %>% left_join(sum,by=c("matrix","row","col")) %>% 
      mutate(par = ifelse(identified,0,ifelse(par==0,par2,par)))
    
  } else {
    # Add free parameter numbers:
    modCopy@parameters$par[modCopy@parameters$par==0 & !modCopy@parameters$identified] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$par==0 & !modCopy@parameters$identified))
    
    # For each group, free all parameters from equality constraints:
    if (type == "free"){
      if (nrow(modCopy@sample@groups)>1){
        for (g in 2:nrow(modCopy@sample@groups)){
          modCopy@parameters$par[modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified))
        }
      }      
    }
    
  }
  
  # Identify:
  # modCopy <- identify(modCopy)
  
  # Remake the model matrix:
  modCopy@extramatrices$M <- Mmatrix(modCopy@parameters)
  
  
  # Check if gradient and hessian are present:
  # gradient <- !is.null(x@fitfunctions$gradient)
  # hessian <- !is.null(x@fitfunctions$hessian)
  # information <- !is.null(x@fitfunctions$information)
  
  # Compute a gradient:
  g <- psychonetrics_gradient(parVector(modCopy), modCopy)
  # if (gradient){
  #   g <- x@fitfunctions$gradient(parVector(modCopy), modCopy)
  # } else {
  #   g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  
  # Compute the Fisher information:
  # if (information){
  #   H <- x@fitfunctions$information(modCopy)
  # } else if (hessian){
  #   H <- x@fitfunctions$hessian(parVector(modCopy), modCopy)
  # } else if (gradient){
  #   H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(modCopy), model=modCopy) 
  # } else {
  #   H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  H <- psychonetrics_FisherInformation(modCopy)
  
  # For every new parameter:
  curMax <- max(x@parameters$par)
  

  # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & x@parameters$par == 0]))){
  # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & (x@parameters$par == 0 | duplicated(x@parameters$par)|rev(duplicated(rev(x@parameters$par))))]))){
  for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par != 0]))){
    ind <- i
    curInds <- seq_len(curMax)
    curInds <- curInds[curInds != i]
    numerator <- as.numeric(H[ind,ind] - H[ind,curInds,drop=FALSE] %*% corpcor::pseudoinverse(H[curInds,curInds]) %*% H[curInds,ind,drop=FALSE])
    
    # Effective n:
    groups <- unique(modCopy@parameters$group_id[modCopy@parameters$par == i])
    Neff <- sum(modCopy@sample@groups$nobs[modCopy@sample@groups$id %in% groups])
    
    
    if (abs(numerator) < 1e-10){
      mi <- NA
      p <- NA
    } else {
      # mi <- (n/Neff) * n * (0.5 * g[i]^2)/numerator
      # mi <- (n/Neff) * ((-n*g[i])^2)/numerator
      mi <- ((-n*g[i])^2)/numerator
      # (n/Neff) * n * (0.5 * (-(n/2)*g[i])^2)/numerator
      
      mi <- as.numeric(mi)
      p <- pchisq(mi,df = 1,lower.tail = FALSE)     
      # epc <- g[i]
    }
    
    
    
    
    if (type == "normal"){
      x@parameters$mi[modCopy@parameters$par == i] <- round(mi,10) # round(mi, 3)
      x@parameters$pmi[modCopy@parameters$par == i] <- round(p,10)
      # x@parameters$epc[modCopy@parameters$par == i] <- round(epc,10)
    } else if (type == "free"){
      x@parameters$mi_free[modCopy@parameters$par == i] <- round(mi,10) # round(mi, 3)
      x@parameters$pmi_free[modCopy@parameters$par == i] <- round(p,10)
    } else {
      x@parameters$mi_equal[modCopy@parameters$par == i] <- round(mi,10) # round(mi,3)
      x@parameters$pmi_equal[modCopy@parameters$par == i] <- round(p, 10)
    }
    
  }
  
  # Make the rest nicer:
  if (type == "normal"){
    x@parameters$mi[is.na(x@parameters$mi)] <- 0
    x@parameters$pmi[is.na(x@parameters$mi)] <- 1
    # x@parameters$epc[is.na(x@parameters$mi)] <- 0
  } else if (type == "free"){
    x@parameters$mi_free[is.na(x@parameters$mi_free)] <- 0
    x@parameters$pmi_free[is.na(x@parameters$mi_free)] <- 1
    # x@parameters$epc_free[is.na(x@parameters$mi_free)] <- 0
  } else {
    x@parameters$mi_equal[is.na(x@parameters$mi_equal)] <- 0
    x@parameters$pmi_equal[is.na(x@parameters$mi_equal)] <- 1
    # x@parameters$epc_equal[is.na(x@parameters$mi_equal)] <- 0
  }
  
  
  return(x)
}