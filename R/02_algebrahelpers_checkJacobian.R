# Simple function to check if a gradient or hessian is correct:
checkJacobian <- function(x, f = psychonetrics_fitfunction, jac = psychonetrics_gradient, transpose = FALSE, plot = TRUE,  perturbStart = FALSE){
  start <- parVector(x)
  
  if (perturbStart){
    start <- start + runif(length(start),0,0.25)
  }
  
  # Analytic:
  analytic <- jac(start, x)
  
  # If not a matrix, make matrix:
  if (!is.matrix(analytic)){
    analytic <- matrix(analytic)
  }
  
  numeric <- numDeriv::jacobian(f,start,model=x)
  
  # If not a matrix, make matrix:
  if (!is.matrix(numeric)){
    numeric <- matrix(numeric)
  }
  
  # transpose:
  if (transpose){
    numeric <- t(numeric)
  }
  
  # plot:
  if (plot){
    plot(Vec(analytic),Vec(numeric),xlab="analytic",ylab="numeric")
    abline(0,1)
  }
  
  return(list(
    analytic = analytic,
    numeric = numeric
  ))
}


# Same, but first replaces observed values with their implied ones
checkFisher <- function(x, f = psychonetrics_gradient, fis = psychonetrics_FisherInformation, transpose = FALSE, plot = TRUE,  perturbStart = FALSE){
  # x <- expectedmodel(x)
  # start <- parVector(x)
  # prep <- prepareModel(start, x)
  # for (g in 1:nrow(x@sample@groups)){
  #   x@sample@means[[g]] <- prep$groupModels[[g]]$mu
  #   x@sample@covs[[g]] <- prep$groupModels[[g]]$sigma
  #   
  #   if (length(x@sample@fimldata) > 0){
  #     nPat <- length(x@sample@fimldata[[g]])
  #     for (i in seq_len(nPat)){
  #       x@sample@fimldata[[g]][[i]]$means <- prep$groupModels[[g]]$mu[!x@sample@fimldata[[g]][[i]]$pattern]
  #       if (!all(x@sample@fimldata[[g]][[i]]$S == 0)){
  #         x@sample@fimldata[[g]][[i]]$S <- prep$groupModels[[g]]$sigma[!x@sample@fimldata[[g]][[i]]$pattern,!x@sample@fimldata[[g]][[i]]$pattern]
  #       }
  #     }
  #   }
  # }
  
  # 
  # 
  # if (perturbStart){
  #   start <- start + runif(length(start),0,0.25)
  # }
  
  # Analytic:
  analytic <- fis(x)
  
  # If not a matrix, make matrix:
  if (!is.matrix(analytic)){
    analytic <- matrix(analytic)
  }
  
  # Numeric:
  # numeric <- 2 * sum(x@sample@groups$nobs) * numDeriv::jacobian(f,start,model=x)
  numeric <- psychonetrics_FisherInformation(x, analytic = FALSE)
  
  # If not a matrix, make matrix:
  if (!is.matrix(numeric)){
    numeric <- matrix(numeric)
  }
  
  # transpose:
  if (transpose){
    numeric <- t(numeric)
  }
  
  # plot:
  if (plot){
    plot(Vec(analytic),Vec(numeric),xlab="analytic",ylab="numeric")
    abline(0,1)
  }
  
  return(list(
    analytic = analytic,
    numeric = numeric
  ))
}