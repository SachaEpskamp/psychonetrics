quantiletransform <- function(x){
  xNoNA <- x[!is.na(x)]
  ord <- order(xNoNA)
  sorted <- sort(x)
  xWithminInf <- c(-Inf,sorted)
  nBelow <- cumsum(xWithminInf[-1] > xWithminInf[-length(xWithminInf)])
  p <- nBelow / (max(nBelow)+1)
  q <- qnorm(p)
  xTrans <- x
  xTrans[!is.na(xTrans)][ord] <- q
  return(xTrans)
}