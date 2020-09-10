# Full function:
addMIs <- function(x, matrices = "all", type =  c("normal","free","equal"),verbose,analyticFisher=TRUE){
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  # full <- TRUE
  # if (full){
  if (verbose){
    message("Computing modification indices...")
  }
  
  tryres <- try({
    
    if ("normal" %in% type){
      
      x <-  x %>% addMIs_inner_full(type = "normal",analyticFisher=analyticFisher)     
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
    
  })
  
  if (is(tryres, "try-error")){
    stop("Failed to compute modification indices. Try a different optimizer with setoptimizer(...) or report the problem on github.com/SachaEpskamp/psychonetrics.")
  }
  
  return(x)
  
}


# Add the modification indices (FULL VERSION):
addMIs_inner_full <- function(x, type =  c("normal","free","equal"),analyticFisher=TRUE){
  type <- match.arg(type)
  
  # If no constrained parameters, nothing to do!
  if (!any(x@parameters$par == 0) & !any(duplicated(x@parameters$par))){
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
      mutate_(par = ~ifelse(identified,0,ifelse(par==0,par2,par)))
    
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
  # modCopy@extramatrices$M <- Mmatrix(modCopy@parameters)
  # 
  if (modCopy@cpp){
    modCopy@extramatrices$M   <- Mmatrix_cpp(modCopy@parameters )
  } else {
    modCopy@extramatrices$M  <- Mmatrix(modCopy@parameters)  
  }
  
  
  # Check if gradient and hessian are present:
  # gradient <- !is.null(x@fitfunctions$gradient)
  # hessian <- !is.null(x@fitfunctions$hessian)
  # information <- !is.null(x@fitfunctions$information)
  
  # Compute a gradient:
  if (modCopy@cpp){
    g <- psychonetrics_gradient_cpp(parVector(modCopy), modCopy)
  } else {
    g <- psychonetrics_gradient(parVector(modCopy), modCopy)
    
  }
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
  # Total sample size:
  nTotal <- sum(x@sample@groups$nobs)
  
  # FIXME: 4 * n could be nicer here probably
  if (modCopy@cpp){
    H <- 4 * nTotal * as(psychonetrics_FisherInformation_cpp(modCopy, analyticFisher), "matrix")
  } else {
    H <- 4 * nTotal * as(psychonetrics_FisherInformation(modCopy, analyticFisher), "matrix")    
  }
  
  
  # For every new parameter:
  curMax <- max(x@parameters$par)
  
  ### NEW
  curInds <- seq_len(curMax)
  newInds <- curMax + seq_len(max(modCopy@parameters$par) - curMax)
  V <-  H[newInds,newInds] - H[newInds,curInds,drop=FALSE] %*% solve_symmetric(H[curInds,curInds]) %*% H[curInds,newInds,drop=FALSE]
  V.diag <- diag(V)
  idx <- which(V.diag < sqrt(.Machine$double.eps))
  if(length(idx) > 0L) {
    V.diag[idx] <- as.numeric(NA)
  }
  # How many in total?
  nTotalPars <- length(c(curInds,newInds))
  
  # 
  # 
  # # Effective N:
  # if (nrow(x@sample@groups) > 1){
  #   par <- modCopy@parameters
  #   # par$id[!modCopy@parameters$identified] <- seq_len(sum(!modCopy@parameters$identified))
  #   par <- par %>% left_join(modCopy@sample@groups, by = c("group_id" = "id")) %>%
  #     group_by(par) %>% summarize(Neff = sum(nobs))
  #   Neff <- numeric(max(par$par))
  #   Neff[par$par[par$par!=0]] <- par$Neff[par$par!=0]
  # } else {
  #   Neff <- rep(x@sample@groups$nobs[1],nTotalPars)
  # }
  
  
  # MIs:
  # All MIs:
  mi <- numeric(nTotalPars)
  
  # mi[newInds] <- ifelse(abs(V.diag) < 1e-10,0,((-Neff[newInds]*g[newInds])^2)/V.diag)
  mi[newInds] <- ifelse(abs(V.diag) < 1e-10,0,((-nTotal*g[newInds])^2)/V.diag)
  if (length(curInds) > 0){
    # mi[curInds] <- ((-Neff[curInds]*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
    mi[curInds] <- ((-nTotal*g[curInds])^2)/diag(H[curInds,curInds,drop=FALSE])
  }
  p <- pchisq(mi,df = 1,lower.tail = FALSE)     
  
  # Compute epc:
  # d <- 0.5 * (-1 * Neff) * g
  d <- 0.5 * (-1 * nTotal) * g
  # needed? probably not; just in case
  d[which(abs(d) < 1e-15)] <- 1.0
  
  # Expected parameter change:
  epc <-   mi / d 
  
  # Which to fill:
  fillInds <- match(c(curInds,newInds),modCopy@parameters$par)
  if (type == "normal"){
    x@parameters$mi[fillInds[!is.na(fillInds)]] <- round(mi[!is.na(fillInds)],10) # round(mi, 3)
    x@parameters$pmi[fillInds[!is.na(fillInds)]] <- round(p[!is.na(fillInds)],10)
    x@parameters$epc[fillInds[!is.na(fillInds)]] <- round(epc[!is.na(fillInds)],10)
  } else if (type == "free"){
    x@parameters$mi_free[fillInds[!is.na(fillInds)]] <- round(mi[!is.na(fillInds)],10) # round(mi, 3)
    x@parameters$pmi_free[fillInds[!is.na(fillInds)]] <- round(p[!is.na(fillInds)],10)
    x@parameters$epc[fillInds[!is.na(fillInds)]] <- round(epc[!is.na(fillInds)],10)
  } else {
    x@parameters$mi_equal[fillInds[!is.na(fillInds)]] <- round(mi[!is.na(fillInds)],10) # round(mi,3)
    x@parameters$pmi_equal[fillInds[!is.na(fillInds)]] <- round(p[!is.na(fillInds)], 10)
    x@parameters$epc[fillInds[!is.na(fillInds)]] <- round(epc[!is.na(fillInds)],10)
  }
  
  return(x)
}