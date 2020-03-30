setoptimizer <- function(x, optimizer = c("default","nlminb","ucminf","cpp_L-BFGS-B","cpp_BFGS","cpp_CG","cpp_SANN","cpp_Nelder-Mead")){
  optimizer <- match.arg(optimizer)
  if (optimizer == "default"){
    optimizer <- defaultoptimizer(x)
  }
  x@optimizer <- optimizer
  
  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@optimizer <- optimizer
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@optimizer <- optimizer
  }
  
  x
}