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
  # } 
  # else {
  #   if (verbose){
  #     message("Computing modification indices...")
  #   }
  #   x <-  x %>% addMIs_inner(equal=FALSE,approxHessian=approxHessian,matrices=matrices) 
  #   if (nrow(x@sample@groups) > 1){
  #     if (verbose){
  #       message("Computing group-constrained indices...")
  #     }
  #     x <- x %>% addMIs_inner(equal=TRUE,approxHessian=approxHessian,matrices=matrices) 
  #   }
  #   return(x)
  # }
}

# 
# # Add the modification indices (ONLY FIXED PARAMETER VERSION):
# addMIs_inner <- function(x, equal = FALSE,approxHessian=FALSE, matrices = "all"){
# 
#   if (any(matrices == "all")){
#     matrices <- x@matrices$name
#   }
#   
#   # If no constrained parameters, nothing to do!
#   if (!any(x@parameters$par == 0)){
#     return(x)
#   }
#   
#   # Clear old MIs:
#   if (!equal){
#     x@parameters$mi[] <- 0
#     x@parameters$pmi[] <- NA
#   } else {
#     x@parameters$mi_equal[] <- 0
#     x@parameters$pmi_equal[] <- NA
#   }
#   
#   # Sample size:
#   n <- sum(x@sample@groups$nobs)
#   
#   # Add two kinds of MIs, one for all fixed parameters free, and one for all fixed free but constrained per group #
#   # Fully free:
#   # Copy the model:
#   modCopy <- x
#   
#   # Obtain the full set of parameters that are constrained across all groups:
#   if (equal){
#     sum <- modCopy@parameters %>% group_by_("matrix","row","col") %>% summarize_(anyConstrained = ~any(fixed)) %>% 
#       filter_(~anyConstrained,~matrix %in% matrices)
#     if (nrow(sum) > 0){
#       # Add a unique number to each:
#       sum$par2 <-  max(modCopy@parameters$par) + seq_len(nrow(sum))
#       
#       # Left join back:
#       modCopy@parameters <- modCopy@parameters %>% left_join(sum,by=c("matrix","row","col")) %>% 
#         mutate(par = ifelse(identified,0,ifelse(par==0,par2,par)))      
#     }
# 
#     
#   } else {
#     # Add free parameter numbers:
#     modCopy@parameters$par[modCopy@parameters$par==0 & !modCopy@parameters$identified & modCopy@parameters$matrix %in% matrices] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$par==0 & !modCopy@parameters$identified& modCopy@parameters$matrix %in% matrices))
#     
#     # For each group, free all parameters from equality constraints (full version only)
#     # if (nrow(modCopy@sample@groups)>1){
#     #   for (g in 2:nrow(modCopy@sample@groups)){
#     #     modCopy@parameters$par[modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified] <- max(modCopy@parameters$par) + seq_len(sum(modCopy@parameters$group_id == g & duplicated(modCopy@parameters$par) & !modCopy@parameters$identified))
#     #   }
#     # }
#   }
# 
#   if (max(modCopy@parameters$par) == max(x@parameters$par)){
#     # Nothing to do
#     return(x)
#   }
#   
#   # Indices of the new parameters:
#   newPars <- max(x@parameters$par) + seq_len(max(modCopy@parameters$par) - max(x@parameters$par))
#   oldPars <- seq_len(max(x@parameters$par))
#   
# 
#   # Remake the model matrix:
#   modCopy@fitfunctions$extramatrices$M <- Mmatrix(modCopy@parameters)
#   
# 
#   
#   # Check if gradient and hessian are present:
#   gradient <- !is.null(x@fitfunctions$gradient)
#   hessian <- !is.null(x@fitfunctions$hessian)
#   
#   if (approxHessian || is.null(x@optim$inverseHessian)){
#     if (hessian){
#       H <- x@fitfunctions$hessian(parVector(x), x)
#     } else if (gradient){
#       H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(x), model=x) 
#     } else {
#       H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(x), model=x) 
#     }
#     
#     # Invert:
#     Hinv <- corpcor::pseudoinverse(H)
#   } else {
#     Hinv <- x@optim$inverseHessian
#   }
#   
#   # Compute a gradient:
#   if (gradient){
#     g <- subGradient(x@fitfunctions$gradient,newPars,x=parVector(modCopy), model =modCopy)
#   } else {
#     g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
#   }
#   
#   
#   # Smart removal of some parameters:
#   # if (all(abs(g) < gradTol)){
#   #   # Nothing to do..
#   #   return(x)
#   # }
#   
#   removePars <- newPars[abs(g) < gradTol]
# 
#   if (length(removePars) > 0){
#     # Update model:
#     modCopy@parameters$par[modCopy@parameters$par %in% removePars] <- 0  
#     modCopy@parameters <- parRelabel(modCopy@parameters)
#     newPars <- max(x@parameters$par) + seq_len(max(modCopy@parameters$par) - max(x@parameters$par))
#     if (gradient){
#       g <- subGradient(x@fitfunctions$gradient,newPars,x=parVector(modCopy), model =modCopy)
#     } else {
#       g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
#     }
#   } 
# 
#   # Compute a hessian:
#   if (hessian){
#     Hnew <- subHessian(hess =  x@fitfunctions$hessian, inds = newPars, x= parVector(modCopy),model = modCopy)
#   } else if (gradient){
#     Hnew <- numDeriv::jacobian(subGradient, parVector(modCopy),gr = x@fitfunctions$gradient, inds = newPars, model=modCopy)
#   } else {
#     Hnew <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy)[newPars,]
#   }
#   
#   # For every new parameter:
#   curMax <- max(x@parameters$par)
#   
#   
#   # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & x@parameters$par == 0]))){
#   # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par>1 & (x@parameters$par == 0 | duplicated(x@parameters$par)|rev(duplicated(rev(x@parameters$par))))]))){
#   # for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par != 0]))){
#   for (i in seq_along(newPars)){
#     mi <- n * (0.5 * g[i]^2)/(Hnew[i,curMax+i] - Hnew[i,oldPars,drop=FALSE] %*% Hinv %*% t(Hnew[i,oldPars,drop=FALSE]))
#     p <- pchisq(mi,df = 1,lower.tail = FALSE)      
#     
#     if (equal){
#       x@parameters$mi_equal[modCopy@parameters$par == newPars[i]] <- round(mi,10) # round(mi,3)
#       x@parameters$pmi_equal[modCopy@parameters$par == newPars[i]] <- round(p, 10)
#     } else {
#       x@parameters$mi[modCopy@parameters$par == newPars[i]] <- round(mi,10) # round(mi, 3)
#       x@parameters$pmi[modCopy@parameters$par == newPars[i]] <- round(p,10)
#     }
#     
#   }
#   
#   if (equal){
#     x@parameters$mi_equal[is.na(x@parameters$mi_equal)] <- 0
#     x@parameters$pmi_equal[is.na(x@parameters$pmi_equal)] <- 1
#   } else {
#     x@parameters$mi[is.na(x@parameters$mi)] <- 0
#     x@parameters$pmi[is.na(x@parameters$pmi)] <- 1
#   }
#   
#   return(x)
# }


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
  modCopy@fitfunctions$extramatrices$M <- Mmatrix(modCopy@parameters)
  
  
  # Check if gradient and hessian are present:
  gradient <- !is.null(x@fitfunctions$gradient)
  hessian <- !is.null(x@fitfunctions$hessian)
  information <- !is.null(x@fitfunctions$information)
  
  # Compute a gradient:
  if (gradient){
    g <- x@fitfunctions$gradient(parVector(modCopy), modCopy)
  } else {
    g <- numDeriv::grad(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  }
  
  # Compute a hessian:
  # if (hessian){
  #   H <- x@fitfunctions$hessian(parVector(modCopy), modCopy)
  # } else if (gradient){
  #   H <- numDeriv::jacobian(x@fitfunctions$gradient,parVector(modCopy), model=modCopy) 
  # } else {
  #   H <- numDeriv::hessian(x@fitfunctions$fitfunction,parVector(modCopy), model=modCopy) 
  # }
  if (information){
    H <- x@fitfunctions$information(modCopy)
  } else if (hessian){
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
  for (i in sort(unique(modCopy@parameters$par[modCopy@parameters$par != 0]))){
    ind <- i
    curInds <- seq_len(curMax)
    curInds <- curInds[curInds != i]
    numerator <- as.numeric(H[ind,ind] - H[ind,curInds,drop=FALSE] %*% solve(H[curInds,curInds]) %*% H[curInds,ind,drop=FALSE])
    
    # Effective n:
    groups <- unique(modCopy@parameters$group_id[modCopy@parameters$par == i])
    Neff <- sum(modCopy@sample@groups$nobs[modCopy@sample@groups$id %in% groups])
    
    
    if (abs(numerator) < 1e-10){
      mi <- NA
      p <- NA
    } else {
      mi <- (n/Neff) * n * (0.5 * g[i]^2)/numerator
      mi <- as.numeric(mi)
      p <- pchisq(mi,df = 1,lower.tail = FALSE)     
      epc <- g[i]
    }
    
    
    
    
    if (type == "normal"){
      x@parameters$mi[modCopy@parameters$par == i] <- round(mi,10) # round(mi, 3)
      x@parameters$pmi[modCopy@parameters$par == i] <- round(p,10)
      x@parameters$epc[modCopy@parameters$par == i] <- round(epc,10)
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