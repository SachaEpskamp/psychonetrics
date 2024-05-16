# Simple function to check if a gradient or hessian is correct:
checkJacobian <- function(x, f = "default", jac = "default", transpose = FALSE, plot = TRUE,  perturbStart = FALSE,method="Richardson"){
  
  start <- parVector(x)
  
  if (identical(f,"default")){
    f <- psychonetrics_fitfunction_cpp
  }
  
  if (identical(jac,"default")){
    jac <- psychonetrics_gradient_cpp
  }

  if (perturbStart){
    start <- start + runif(length(start),0,0.25)
  }
  
  

  # Analytic:
  analytic <- jac(start, x)
  
  
  # If not a matrix, make matrix:
  if (!is.matrix(analytic)){
    analytic <- matrix(analytic)
  }

  numeric <- numDeriv::jacobian(f,start,model=x,method=method,method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))

  
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
    plot(as.vector(as.numeric(Vec(analytic))),as.vector(as.numeric(Vec(numeric))),xlab="analytic",ylab="numeric")
    abline(0,1)
  }
  

  return(list(
    analytic = as.vector(analytic),
    numeric = as.vector(numeric)
  ))
}


# Same, but first replaces observed values with their implied ones
checkFisher <- function(x, f = "default", fis = "default", transpose = FALSE, plot = TRUE,  perturbStart = FALSE){
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
  
  if (identical(f,"default")){
    f <- psychonetrics_fitfunction_cpp
  }
  
  if (identical(fis,"default")){
    fis <- psychonetrics_FisherInformation_cpp
  }
  
  
  
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
    analytic = as.matrix(analytic),
    numeric = as.matrix(numeric)
  ))
}