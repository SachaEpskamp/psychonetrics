# Run a model!
runmodel <- function(
  x, # psychonetrics model
  # stepwise = FALSE, # Stepwise up search with modification indices?
  matrices, # Matrices to search
  level = c("gradient","fitfunction"),
  addfit = TRUE,
  addMIs = TRUE,
  addSEs=TRUE,
  log = TRUE,
  verbose = TRUE,
  # optimizer = c("default","ucminf","nlminb"),
  optim.control = list()
  # inverseHessian = TRUE
){
  optimizer <- x@optimizer
  # Default:
  if (optimizer == "default"){
    if (any(grepl("omega",x@parameters$matrix))){
      optimizer <- "nlminb"      
    } else {
      optimizer <- "ucminf"
    }

  }
  
  
  # optimizer <- match.arg(optimizer)
  level <- match.arg(level)
  if (!is(x,"psychonetrics")){
    stop("input is not a 'psychonetrics' object")
  }
  
  
  # Evaluate baseline model:
  if (!is.null(x@baseline_saturated$baseline) && !x@baseline_saturated$baseline@computed){
    if (verbose) message("Estimating baseline model...")
    # Run:
    
    x@baseline_saturated$baseline <- runmodel(x@baseline_saturated$baseline, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE)
  }
  
  # Evaluate saturated model:
  
  if (!is.null(x@baseline_saturated$saturated) && !x@baseline_saturated$saturated@computed){
    if (verbose) message("Estimating saturated model...")
    # Run:
    x@baseline_saturated$saturated <- runmodel(x@baseline_saturated$saturated, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE)
  }
  
  
  
  # # nlminb control pars:
  # control.nlminb <- list(eval.max=20000L,
  #                        iter.max=10000L,
  #                        trace=0L,
  #                        #abs.tol=1e-20, ### important!! fx never negative
  #                        abs.tol=(.Machine$double.eps * 10),
  #                        # rel.tol=1e-10,
  #                        rel.tol=1e-5,
  #                        #step.min=2.2e-14, # in =< 0.5-12
  #                        step.min=1.0, # 1.0 in < 0.5-21
  #                        step.max=1.0,
  #                        x.tol=1.5e-8,
  #                        xf.tol=2.2e-14)
  # control.nlminb <- modifyList(control.nlminb, nlminb.control)
  # control <- control.nlminb[c("eval.max", "iter.max", "trace",
  #                             "step.min", "step.max",
  #                             "abs.tol", "rel.tol", "x.tol", "xf.tol")]
  
  # Start and bounts:
  start <- parVector(x)
  lower <- lowerBound(x)
  upper <- upperBound(x)
  
  # Check if Gradient and hessian are present:
  # if (level == "default"){
  #   if (!is.null(x@fitfunctions$gradient) & !is.null(x@fitfunctions$hessian)){
  #     level <- "hessian"
  #   } else  if (!is.null(x@fitfunctions$gradient)){
  #     level <- "gradient"
  #   } else {
  #     level <- "fitfunction"
  #   }    
  # }
  
  # Default optimizer:
  # optimizer <- match.arg(optimizer)
  # if (optimizer == "default"){
  #   if (level == "hessian"){
  #     optimizer <- "nlminb"
  #   } else {
  #     optimizer <- "ucminf"
  #   }
  #   
  # }
  
  # if (optimizer%in% c("Nelder-Mead","L-BFGS-B","ucminf") & level == "hessian"){
  # if (optimizer%in% c("ucminf") & level == "hessian"){
  #   warning("Optimizer does not support analytical Hessian. Using numeric Hessian instead.")
  #   level <- "gradient"
  # }
  
  if (verbose) message("Estimating model...")
  # Form optimizer arguments:
  
  # if (optimizer == "nlminb"){
  #   control.nlminb <- list(eval.max=20000L,
  #                          iter.max=10000L,
  #                          trace=0L,
  #                          #abs.tol=1e-20, ### important!! fx never negative
  #                          abs.tol=(.Machine$double.eps * 10),
  #                          # rel.tol=1e-10,
  #                          rel.tol=1e-5,
  #                          #step.min=2.2e-14, # in =< 0.5-12
  #                          step.min=1.0, # 1.0 in < 0.5-21
  #                          step.max=1.0,
  #                          x.tol=1.5e-8,
  #                          xf.tol=2.2e-14)
  # 
  #   optim.control$control <- control.nlminb
  #   optim.control$start <- start
  #   optim.control$objective <- psychonetrics_fitfunction
  #   optim.control$lower <- lower
  #   optim.control$upper <- upper
  #   optim.control$model <- x
  #   if (level != "fitfunction"){
  #     # optim.control$gradient <- x@fitfunctions$gradient
  #     optim.control$gradient <- psychonetrics_gradient
  #   }
  #   # if (level == "hessian"){
  #   #   x@fitfunctions$hessian <- x@fitfunctions$hessian
  #   # }
  #   
  #   # Run model:
  #   optim.out <- do.call(nlminb,optim.control)
  #     
  #   # Update model:
  #   x <- updateModel(optim.out$par,x)
  #   
  #   # Make list:
  #   optimresults <- list(
  #     par = optim.out$par,
  #     value = optim.out$objective,
  #     # iterations = optim.out$iterations,
  #     message = optim.out$message,
  #     optimizer = optimizer
  #   )
  #   
  #   # # Compute inverse Hessian if needed:
  #   # if (inverseHessian){
  #   #   if (level == "hessian"){
  #   #     H <- x@fitfunctions$hessian(optim.out$par, x)
  #   #   } else if (level == "gradient"){
  #   #     H <- numDeriv::jacobian(x@fitfunctions$gradient, optim.out$par, model = x)
  #   #   } else {
  #   #     H <- numDeriv::hessian(x@fitfunctions$fitfunction, optim.out$par, model = x)        
  #   #   }
  #   #   Hinv <- corpcor::pseudoinverse(H)
  #   #   optimresults$inverseHessian <- Hinv
  #   # }
  # } else if (optimizer == "ucminf"){
  #   optim.control$par <- start
  #   optim.control$fn <- psychonetrics_fitfunction
  #   optim.control$model <- x
  #   if (level != "fitfunction"){
  #     optim.control$gr <- psychonetrics_gradient
  #   }
  #   # if (inverseHessian){
  #   #   optim.control$hessian <- 2
  #   # }
  #   
  #   # Run model:
  #   optim.out <- do.call(ucminf,optim.control)
  #   
  #   # Update model:
  #   x <- updateModel(optim.out$par,x)
  #   
  #   # Make list:
  #   optimresults <- list(
  #     par = optim.out$par,
  #     value = optim.out$value,
  #     message = optim.out$message,
  #     optimizer = optimizer
  #   )
  #   
  #   # if (inverseHessian){
  #   #   optimresults$inverseHessian <- optim.out$invhessian
  #   # }
  # }
  
  
  optim.control$par <- start
  optim.control$fn <- psychonetrics_fitfunction
  optim.control$model <- x
  if (level != "fitfunction"){
    optim.control$gr <- psychonetrics_gradient
  }
  
  # Add method:
  optim.control$method <- optimizer
  
  # Add bounds:
  if (optimizer %in% c("nlmimb","L-BFGS-B","lbfgs")){
    optim.control$lower <- lower
    optim.control$upper <- upper
  }
  
  # Run model:
  optim.out <- do.call(optimr,optim.control)
  
  # Update model:
  x <- updateModel(optim.out$par,x,updateMatrices = TRUE)
  
  # Make list:
  # optimresults <- list(
  #   par = optim.out$par,
  #   value = optim.out$value,
  #   message = optim.out$message,
  #   optimizer = optimizer
  # )
  optimresults <- optim.out
  optimresults$optimizer <- optimizer
  
  # 
  #   ### START OPTIMIZATION ###
  #   if (level == "fitfunction"){
  #     if (optimizer == "nlminb"){
  #       optim.out <- nlminb(start=start,
  #                           objective=x@fitfunctions$fitfunction,
  #                           gradient=NULL,
  #                           hessian = NULL,
  #                           lower=lower,
  #                           upper=upper,
  #                           model = x,
  #                           control=control
  #       )      
  #     } else if (optimizer == "ucminf"){
  #       out <- ucminf(start, x@fitfunctions$fitfunction, model = x)
  #       optim.out <- list(
  #         par = out$par,
  #         objective = out$value,
  #         convergence = out$convergence,
  #         message = out$message
  #       )
  #     } else {
  #       out <- optim(start,
  #                    fn=x@fitfunctions$fitfunction,
  #                    model=x,
  #                    lower=lower,
  #                    upper = upper,
  #                    # gr = x@fitfunctions$gradient,
  #                    method = optimizer)
  #       optim.out <- list(
  #         par = out$par,
  #         objective = out$value,
  #         convergence = out$convergence,
  #         message = out$message
  #       )
  #     }
  #     
  #     # scale=SCALE, # FIXME: What is this in lavaan?
  #   } else if (level == "gradient"){
  #     if (optimizer == "nlminb"){
  #       optim.out <- nlminb(start=start,
  #                           objective=x@fitfunctions$fitfunction,
  #                           gradient=x@fitfunctions$gradient,
  #                           hessian = NULL,
  #                           lower=lower,
  #                           upper=upper,
  #                           model = x,
  #                           control=control
  #       )
  #     } else if (optimizer == "ucminf"){
  #       out <- ucminf(start, x@fitfunctions$fitfunction,gr = x@fitfunctions$gradient, model = x)
  #       optim.out <- list(
  #         par = out$par,
  #         objective = out$value,
  #         convergence = out$convergence,
  #         message = out$message
  #       )
  #     } else {
  #       out <- optim(start,
  #                    fn=x@fitfunctions$fitfunction,
  #                    model=x,
  #                    lower=lower,
  #                    upper = upper,
  #                    gr = x@fitfunctions$gradient,
  #                    method = optimizer)
  #       optim.out <- list(
  #         par = out$par,
  #         objective = out$value,
  #         convergence = out$convergence,
  #         message = out$message
  #       )
  #     }
  #     
  #   } else {
  #     
  #     optim.out <- nlminb(start=start,
  #                         objective=x@fitfunctions$fitfunction,
  #                         gradient=x@fitfunctions$gradient,
  #                         hessian = x@fitfunctions$hessian,
  #                         lower=lower,
  #                         upper=upper,
  #                         model = x,
  #                         control=control
  #     )
  #     
  #   }
  #   
  
  x@optim <- optimresults
  x@computed <- TRUE
  x@objective <- optimresults$value
  
  # Add information:
  # if (!is.null(x@fitfunctions$information)){
  x@information <- psychonetrics_FisherInformation(x)
  
  if (!all(x@information == t(x@information))){
    x@information <- 0.5 * (x@information + t(x@information))
  }
  if (any(Re(eigen(x@information)$values) < 0)){
    warning("Information matrix is not positive semi-definite. Model might not be identified.")
  }
  # }
  # Add fit:
  if (addfit){
    x <- addfit(x)
  }
  # Add MIs:
  if (addMIs){
    x <- addMIs(x) 
  }
  # Add SEs:
  if (addSEs){
    x <- addSEs(x)
  }
  
  if (log){
    # Add log:
    x <- addLog(x, "Evaluated model")    
  }
  
  
  # Return model:
  return(x)
}