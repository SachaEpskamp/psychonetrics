usecpp <- function(x, use = TRUE){
  x@cpp <- use
  
  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@cpp <- use
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@cpp <- use
  }
  
  x
}