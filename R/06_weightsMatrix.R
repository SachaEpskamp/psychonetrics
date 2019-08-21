LS_weightsmat <- function(dat, type = c("full","diagonal"), meanstructure = TRUE){
  type <- match.arg(type)
  nvar <- ncol(dat)
  ncase <- nrow(dat)
  
  if (type == "full"){
    Wmat <- WLS_wmat(
      as.matrix(dat),
      colMeans(dat),
      ncase,
      nvar)    
  } else if (type == "diagonal"){
    Wmat <- DWLS_wmat(
      as.matrix(dat),
      colMeans(dat),
      ncase,
      nvar)  
  }

  # If the mean structure is ignored, set the mean part to arbitrary dummy identity matrix:
  if (!meanstructure){
    Wmat[seq_len(nvar),] <- 0
    Wmat[,seq_len(nvar)] <- 0
    Wmat[seq_len(nvar),seq_len(nvar)] <- diag(nvar)
  }
  
  WmatInv <- as(solve_symmetric(as(Wmat,"Matrix")),"Matrix")
  WmatInv
}