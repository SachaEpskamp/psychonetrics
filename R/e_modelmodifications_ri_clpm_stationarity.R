
# FUnction to make the RI CLPM model more stationary:
ri_clpm_stationary <- function(
    x,
    stationary = c( "intercepts", "contemporaneous", "innovation", "temporal")
){
  stopifnot(is(x, "psychonetrics"))
  
  # Check stationary 
  stationary <- match.arg(stationary)
  
  # Verbose:
  verbose <- x@verbose
  
  # Extra matrices:
  vars <- x@extramatrices$vars
  varsDF <- x@extramatrices$varsDF
  
  # Number of variables:
  nVar <- nrow(vars)
  
  # Number of measurements:
  nTime <- ncol(vars)
  
  # Number of total observed:
  n_observed <- nVar * nTime
  
  # Set stationary:
  if (stationary == "intercepts"){
    
    # Average intercept:
    avg_int <- rowMeans(matrix(getmatrix(x, "nu", group = 1),nrow=nVar))
    
    # Fix all nu parameters:
    x <- fixpar(x, matrix = "nu", row = seq_len(n_observed), col = 1)
    
    # Free nu_eta parameters for random intercepts:
    x <- freepar(x, matrix = "nu_eta", row = n_observed + seq_len(nVar), col = 1)
    
    # Set starting values:
    x@parameters$est[x@parameters$matrix == "nu_eta" & !x@parameters$fixed] <- avg_int
    
  } else if (stationary == "temporal") {
    
    # Beta parameter numbers:
    parNums <- which(x@parameters$matrix == "beta" & !x@parameters$fixed)
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar^2){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime-1)
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  } else if (stationary == "contemporaneous") {

    # Beta parameter numbers:
    parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row != x@parameters$col))
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar*(nVar-1)/2){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime+1)
    
    # cut out time = 1 and RI:
    parNums <- parNums[,-c(1,ncol(parNums)),drop=FALSE]
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  }  else if (stationary == "innovation") {
    
    # Beta parameter numbers:
    parNums <- which(grepl("_zeta", x@parameters$matrix) & !x@parameters$fixed & (x@parameters$row == x@parameters$col))
    
    # check if not pruned:
    if (length(parNums) < (nTime-1) * nVar){
      stop("Stationary temporal effects for pruned models are not yet supported")
    }
    
    # collect in matrix:
    parNums <- matrix(parNums,ncol=nTime+1)
    
    # cut out time = 1 and RI:
    parNums <- parNums[,-c(1,ncol(parNums)),drop=FALSE]
    
    # Constrain equal:
    for (i in seq_len(nrow(parNums))){
      x <- parequal(x, parNums[i,])
    }
    
  } else stop(paste0("stationary = '",stationary,"' not implemented!"))
  
  
  # Return model:
  return(x)
}