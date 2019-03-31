# This function will make a dummy multiple group model for all possible missingness patterns
missingpatterns <- function(dat, verbose = TRUE){
  if (verbose){
    message("Computing missingness patterns...")
  }
  
  # Remove rows with full missing:
  dat <- dat[rowSums(is.na(dat)) < ncol(dat),]
  
  
  # Create a dummy dataset with only missings:
  mis <- as(is.na(dat),"sparseMatrix")
  
  # Unique patterns:
  unMis <- mgcv::uniquecombs(mis)

  if (is(unMis,"logical")){
    unMis <- t(as.matrix(unMis))
  }

  # DUmmy sigma for indices:
  nvar <- ncol(dat)
  dumSig <- matrix(0,nvar,nvar)
  dumSig[lower.tri(dumSig,diag=TRUE)] <- nvar + seq_len(sum(lower.tri(dumSig,diag=TRUE)))
  patterns <- vector("list",nrow(unMis))
  
  if (verbose){
    pb <- txtProgressBar(min=0,max=nrow(unMis),style = 3)
  }
  # for every pattern:
  for (i in 1:nrow(unMis)){
    patterns[[i]]$inds <- which(colSums(t(mis) == unMis[i,]) == ncol(mis))
    patterns[[i]]$n <- length(patterns[[i]]$inds)
    patterns[[i]]$pattern <- unMis[i,]
    subDat <- as.matrix(dat[patterns[[i]]$inds,patterns[[i]]$pattern!=1,drop=FALSE])
    patterns[[i]]$means <- colMeans(subDat)
    patterns[[i]]$S <- 1 / patterns[[i]]$n  *  t(subDat) %*% subDat -patterns[[i]]$means %*% t(patterns[[i]]$means)
    
    # FIXME: Obs vec for RcppArma:
    patterns[[i]]$obs <- as.vector(!unMis[i,])
    
    # Means elimination matrix:
    obs <- !patterns[[i]]$pattern
    
    # Indices:
    inds <- c(
      which(obs), # Mean part
      c(dumSig[obs,obs,drop=FALSE]))
    inds <- inds[inds!=0]
    
    # Elimintation matrix:
    patterns[[i]]$L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),nvar + nvar*(nvar+1)/2))

    
    # Duplication matrix: 
    patterns[[i]]$D <- duplicationMatrix(sum(obs))
    
    
    # Stuff that Armadillo understands:
    patterns[[i]]$L <- as( patterns[[i]]$L , "dgCMatrix")
    patterns[[i]]$D <- as( patterns[[i]]$D , "dgCMatrix")
    
    # patterns[[i]]$Lmu <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(dat)))
    
    # Find the proper elimination matrix:
    # inds <- c(dumSig[obs,obs,drop=FALSE])
    # patterns[[i]]$Lsig <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(dumSig)^2))
    if (verbose){
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose){
    close(pb)
  }
  
  patterns
}


