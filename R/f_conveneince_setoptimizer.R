setoptimizer <- function(x, optimizer = c("default","nlminb","ucminf","nloptr_TNEWTON","LBFGS++"),
                         optim.args){
  optimizer <- match.arg(optimizer)
  if (optimizer == "default"){
    optimizer <- defaultoptimizer(x)
  }
  x@optimizer <- optimizer
  
  if (!missing(optim.args)){
    x@optim.args <- optim.args
  }
  
  if (!is.null(x@baseline_saturated$baseline)){
    x@baseline_saturated$baseline@optimizer <- optimizer
    if (!missing(optim.args)){
      x@baseline_saturated$baseline@optim.args <- optim.args
    }
    
  }
  if (!is.null(x@baseline_saturated$saturated)){
    x@baseline_saturated$saturated@optimizer <- optimizer
    
    if (!missing(optim.args)){
      x@baseline_saturated$saturated@optim.args <- optim.args
    }
  }
  
  x
}