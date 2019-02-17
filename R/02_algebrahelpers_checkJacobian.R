# Simple function to check if a gradient or hessian is correct:
checkJacobian <- function(x, f, jac, transpose = FALSE, plot = TRUE,  perturbStart = FALSE){
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
  
  # Numeric:
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