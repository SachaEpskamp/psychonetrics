LS_weightsmat <- function(dat, type = c("full","diagonal")){
  type <- match.arg(type)
  nvar <- ncol(dat)
  ncase <- nrow(dat)
  
  WLS_wmat(
     as.matrix(dat),
   colMeans(dat),
    ncase,
    nvar)
}