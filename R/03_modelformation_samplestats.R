# Inner function to make sample stats object:
samplestats <- function(
  ..., rawts = FALSE){
  if (rawts){
    return(samplestats_rawts(...))
  } else {
    return(samplestats_norawts(...))
  }
}