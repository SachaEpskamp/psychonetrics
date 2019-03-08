setoptimizer <- function(x, optimizer){
  x@optimizer <- optimizer
  
  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@optimizer <- optimizer
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@optimizer <- optimizer
  }
  
  x
}