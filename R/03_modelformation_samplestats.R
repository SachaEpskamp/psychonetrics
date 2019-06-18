# Inner function to make sample stats object:
samplestats <- function(
  ..., rawts = FALSE){
  if (rawts){
    stop("rawts currently not supported")
  } else {
    return(samplestats_norawts(...))
  }
}