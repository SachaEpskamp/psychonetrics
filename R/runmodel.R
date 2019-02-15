# Run a model!
runmodel <- function(
  x, # psychonetrics model
  stepwise = FALSE, # Stepwise up search with modification indices?
  alpha = 0.01, # Alpha to use for modification indices
  matrices, # Matrices to search
  nlminb.control = list(),
  level = c("default","fitfunction","gradient","hessian")
){
  level <- match.arg(level)
  if (!is(x,"psychonetrics")){
    stop("input is not a 'psychonetrics' object")
  }
  if (stepwise){
    stop("'stepwise' is not yet implemented")
  }
  
  # nlminb control pars:
  control.nlminb <- list(eval.max=20000L,
                         iter.max=10000L,
                         trace=0L,
                         #abs.tol=1e-20, ### important!! fx never negative
                         abs.tol=(.Machine$double.eps * 10),
                         rel.tol=1e-10,
                         #step.min=2.2e-14, # in =< 0.5-12
                         step.min=1.0, # 1.0 in < 0.5-21
                         step.max=1.0,
                         x.tol=1.5e-8,
                         xf.tol=2.2e-14)
  control.nlminb <- modifyList(control.nlminb, nlminb.control)
  control <- control.nlminb[c("eval.max", "iter.max", "trace",
                              "step.min", "step.max",
                              "abs.tol", "rel.tol", "x.tol", "xf.tol")]
  
  start <- parVector(x)
  
  # Check if Gradient and hessian are present:
  if (level == "default"){
    if (!is.null(x@fitfunctions$gradient) & !is.null(x@fitfunctions$hessian)){
      level <- "hessian"
    } else  if (!is.null(x@fitfunctions$gradient)){
      level <- "gradient"
    } else {
      level <- "fitfunction"
    }    
  }
  
  ### START OPTIMIZATION ###
  if (level == "fitfunction"){
    optim.out <- nlminb(start=start,
                        objective=x@fitfunctions$fitfunction,
                        gradient=NULL,
                        hessian = NULL,
                        lower=-Inf,
                        upper=Inf,
                        model = x,
                        control=control)
    # scale=SCALE, # FIXME: What is this in lavaan?
  } else if (level == "gradient"){
    optim.out <- nlminb(start=start,
                        objective=x@fitfunctions$fitfunction,
                        gradient=x@fitfunctions$gradient,
                        hessian = NULL,
                        lower=-Inf,
                        upper=Inf,
                        model = x,
                        control=control)
  } else {
    optim.out <- nlminb(start=start,
                        objective=x@fitfunctions$fitfunction,
                        gradient=x@fitfunctions$gradient,
                        hessian = x@fitfunctions$hessian,
                        lower=-Inf,
                        upper=Inf,
                        model = x,
                        control=control)
  }
  
  # Update model:
  x <- updateModel(optim.out$par,x)
  x@optim <- optim.out
  x@computed <- TRUE
  x@objective <- optim.out$objective

  # Return model:
  return(x)
}