generate <- function(x, n = 500){
  if (!is(x,"psychonetrics")){
    stop("Input must be a psychonetrics object")
  }
  
  if (x@distribution != "Gaussian"){
    stop("Only Gaussian data supported at the moment.")
  }
  
  if (x@model == "var1"){
    stop("VAR models not supported yet.")
  }
  
  generate_gaussian(x,n)
}

generate_gaussian <- function(x,n){
  nGroup <- length(x@sample@groups$id)
  data <- do.call(rbind,lapply(seq_len(nGroup),function(i){
    data <- rmvnorm(n, as.vector(x@modelmatrices[[i]]$mu),
                    as.matrix(x@modelmatrices[[i]]$sigma))
    data <- as.data.frame(data)
    names(data) <- x@sample@variables$label
    data$group <- i
    data
  }))
  if (nGroup == 1){
    data <- data[,-ncol(data)]
  }  
  return(data)
}
