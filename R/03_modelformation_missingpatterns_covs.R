# This function will make a dummy multiple group model for all possible missingness patterns
missingpatterns_covs <- function(means, covs, nobs, verbose = TRUE){
  if (verbose){
    message("Computing missingness patterns...")
  }
  # browser()
  # # Remove rows with full missing:
  # dat <- dat[rowSums(is.na(dat)) < ncol(dat),]
  # 
  # 
  # # Create a dummy dataset with only missings:
  # mis <- as(is.na(dat),"sparseMatrix")
  # 
  # # Unique patterns:
  # unMis <- mgcv::uniquecombs(mis)
  # 
  # if (is(unMis,"logical")){
  #   unMis <- t(as.matrix(unMis))
  # }
  # 
  # # DUmmy sigma for indices:
  nvar <- ncol(covs)
  dumSig <- matrix(0,nvar,nvar)
  dumSig[lower.tri(dumSig,diag=TRUE)] <- nvar + seq_len(sum(lower.tri(dumSig,diag=TRUE)))
  # patterns <- vector("list",nrow(unMis))
  # 
  # if (verbose){
  #   pb <- txtProgressBar(min=0,max=nrow(unMis),style = 3)
  # }
  # for every pattern:
  
  
  # Dummy, there is only one pattern!
  patterns <- list(list())
  
  # Obtain the pattern
  misPattern <- rowSums(is.na(covs)) == nrow(covs) & colSums(is.na(covs)) == ncol(covs)
  
  # for (i in 1:nrow(unMis)){
    patterns[[1]]$inds <- seq_len(nobs)
    patterns[[1]]$n <- nobs
    patterns[[1]]$pattern <- misPattern
    patterns[[1]]$means <- as.vector(means[!misPattern])
    patterns[[1]]$S <- as.matrix(covs[!misPattern,!misPattern])
    
    # FIXME: Obs vec for RcppArma:
    patterns[[1]]$obs <- as.vector(!misPattern)
    
    # Means elimination matrix:
    obs <- !patterns[[1]]$pattern
    
    # Indices:
    inds <- c(
      which(obs), # Mean part
      c(dumSig[obs,obs,drop=FALSE]))
    inds <- inds[inds!=0]
    
    # Elimintation matrix:
    patterns[[1]]$L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),nvar + nvar*(nvar+1)/2))

    
    # Duplication matrix: 
    patterns[[1]]$D <- duplicationMatrix(sum(obs))
    
    
    # Stuff that Armadillo understands:
    patterns[[1]]$L <- as( patterns[[1]]$L , "dgCMatrix")
    patterns[[1]]$D <- as( patterns[[1]]$D , "dgCMatrix")
    
    # patterns[[1]]$Lmu <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(dat)))
    
    # Find the proper elimination matrix:
    # inds <- c(dumSig[obs,obs,drop=FALSE])
    # patterns[[1]]$Lsig <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),ncol(dumSig)^2))
    # if (verbose){
    #   setTxtProgressBar(pb, i)
    # }
  # }
  
  # if (verbose){
  #   close(pb)
  # }
  
  patterns
}


