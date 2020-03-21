setoptimizer <- function(x, optimizer = c("cpp_CG","cpp_BFGS","cpp_L-BFGS-B","cpp_SANN","cpp_Nelder-Mead","nlminb","ucminf")){
  optimizer <- match.arg(optimizer)
  x@optimizer <- optimizer
  
  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@optimizer <- optimizer
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@optimizer <- optimizer
  }
  
  x
}