quantiletransform <- function(x){
  xNoNA <- x[!is.na(x)]
  ord <- order(xNoNA)
  sorted <- sort(xNoNA)
  nBelow <- rank(sorted, ties.method = "min")
  p <- nBelow / (max(nBelow)+1)
  q <- qnorm(p)
  xTrans <- x
  xTrans[!is.na(xTrans)][ord] <- q
  return(xTrans)
}