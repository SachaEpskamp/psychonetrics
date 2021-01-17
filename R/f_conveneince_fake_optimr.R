# Fake optimr function that is needed for nlminb because of bug in optimr source code
optimr_fake <- function(...){
  dots <- list(...)
  if (dots$method != "nlminb"){
    return(optimr(...))
  } else {
    
    newargs <- dots
    
    rename_element <- function(x, from, to){
      x[[to]] <- x[[from]]
      x[[from]] <- NULL
      return(x)
    }
    newargs$method <- NULL
    newargs <- rename_element(newargs,"gr","gradient")
    newargs <- rename_element(newargs,"hess","hessian")
    newargs <- rename_element(newargs,"par","start")
    newargs <- rename_element(newargs,"fn","objective")
    
    # Run nlminb:
    res <- do.call(nlminb, newargs)
    
    # Now rename some results:
    res <- rename_element(res,"objective", "value")
    res <- rename_element(res,"evaluations", "counts")
    return(res)
  }
  
  
}