LS_weightsmat <- function(dat, type = c("full","diagonal")){
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

  
  WmatInv <- as(solve_symmetric(as(Wmat,"Matrix")),"Matrix")
  WmatInv
}