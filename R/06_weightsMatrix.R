LS_weightsmat <- function(dat, type = c("full","diagonal"), meanstructure = TRUE, corinput = FALSE){
  type <- match.arg(type)
  nvar <- ncol(dat)
  ncase <- nrow(dat)
  
  if (type == "full"){
    Wmat <- WLS_wmat(
      as.matrix(dat),
      colMeans(dat, na.rm = TRUE),
      ncase,
      nvar)    
  } else if (type == "diagonal"){
    Wmat <- DWLS_wmat(
      as.matrix(dat),
      colMeans(dat, na.rm = TRUE),
      ncase,
      nvar)  
  }

  # If the mean structure is ignored, remove from ACOV matrix
  # FIXME: Nicer to never compute this in the first place!
  if (!meanstructure){
    Wmat <- Wmat[-seq_len(nvar),-seq_len(nvar)]
    
    # Wmat[seq_len(nvar),] <- 0
    # Wmat[,seq_len(nvar)] <- 0
    # Wmat[seq_len(nvar),seq_len(nvar)] <- diag(nvar)
  }
  
  # If corinput, remove variances from the Wmat. Note: only happens when data are standardized.
  # FIXME: Nicer to never compute this in the first place!
  if (corinput){
    inds <- meanstructure * nvar + which(diag(nvar)[lower.tri(diag(nvar),diag=TRUE)] == 1)
    Wmat <- Wmat[-inds,-inds]
  }
  
  WmatInv <- as(solve_symmetric(as(Wmat,"matrix")),"matrix")
  WmatInv
}