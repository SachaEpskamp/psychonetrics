# This function will make a dummy multiple group model for all possible missingness patterns
missingpatterns <- function(dat, verbose = TRUE){
  if (verbose){
    message("Computing missingness patterns...")
  }
  
  # Remove rows with full missing:
  dat <- dat[rowSums(is.na(dat)) < ncol(dat),]
  
  
  # Create a dummy dataset with only missings:
  mis <- as(is.na(dat),"matrix")
  
  # Unique patterns:
  unMis <- mgcv::uniquecombs(mis)

  if (is(unMis,"logical")){
    unMis <- t(as.matrix(unMis))
  }

  # DUmmy sigma for indices:
  nvar <- ncol(dat)
  dumSig <- matrix(0,nvar,nvar)
  dumSig[lower.tri(dumSig,diag=TRUE)] <- nvar + seq_len(sum(lower.tri(dumSig,diag=TRUE)))
  browser
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
    patterns[[i]]$L <- as(patterns[[i]]$L, "dMatrix")
    
    # Duplication matrix: 
    patterns[[i]]$D <- duplicationMatrix(sum(obs))
    # patterns[[i]]$D <- as(as.matrix(patterns[[i]]$D), "dMatrix")
    
    # Stuff that Armadillo understands:
    # Not needed, Arma already understands!
    # patterns[[i]]$L <- as( patterns[[i]]$L , "dMatrix")
    # patterns[[i]]$D <- as( patterns[[i]]$D , "dMatrix")
    
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


# FIXME: Function below is just a copy of the one above and can be much faster ....
# Function for fimlData of all individual rows:
fullfimldata <- function(dat, verbose = TRUE){
  if (verbose){
    message("Storing FIML observations...")
  }
  
  # Remove rows with full missing:
  dat <- dat[rowSums(is.na(dat)) < ncol(dat),]
  
  
  # Create a dummy dataset with only missings:
  mis <- as(is.na(dat),"dMatrix")
  
  # DUmmy sigma for indices:
  nvar <- ncol(dat)
  dumSig <- matrix(0,nvar,nvar)
  dumSig[lower.tri(dumSig,diag=TRUE)] <- nvar + seq_len(sum(lower.tri(dumSig,diag=TRUE)))
  
  patterns <- vector("list",nrow(mis))
  
  if (verbose){
    pb <- txtProgressBar(min=0,max=nrow(mis),style = 3)
  }
  # for every pattern:
  for (i in 1:nrow(mis)){
    patterns[[i]]$inds <- i
    patterns[[i]]$n <- 1 # FIXME: silly to treat every observation as a n=1 group
    patterns[[i]]$pattern <- mis[i,]
    subDat <- as.matrix(dat[patterns[[i]]$inds,patterns[[i]]$pattern!=1,drop=FALSE])
    patterns[[i]]$means <- colMeans(subDat)
    patterns[[i]]$S <- 1 / patterns[[i]]$n  *  t(subDat) %*% subDat - patterns[[i]]$means %*% t(patterns[[i]]$means)
    
    # FIXME: Obs vec for RcppArma:
    patterns[[i]]$obs <- as.vector(!mis[i,])
    
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
    patterns[[i]]$L <- as( patterns[[i]]$L , "dMatrix")
    # patterns[[i]]$D <- as( patterns[[i]]$D , "dMatrix")
    
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


