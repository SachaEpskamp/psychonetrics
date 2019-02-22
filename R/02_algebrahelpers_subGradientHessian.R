# Silly functions to obtain only a subset of the gradient or hessian:
subGradient <- function(gr,inds,...){
  gr(...)[inds]
}

subHessian <- function(hess,inds,...){
  hess(...)[inds,,drop=FALSE]
}