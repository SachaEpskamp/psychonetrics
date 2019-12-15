# Run a model!
runmodel <- function(
  x, # psychonetrics model
  # stepwise = FALSE, # Stepwise up search with modification indices?
  level = c("gradient","fitfunction"),
  addfit = TRUE,
  addMIs = TRUE,
  addSEs=TRUE,
  addInformation = TRUE,
  log = TRUE,
  verbose = TRUE,
  # optimizer = c("default","ucminf","nlminb"),
  optim.control = list(),
  maxtry = 5,
  analyticFisher = TRUE
  # inverseHessian = TRUE
){
  # first check if there are any free parameters:
  if (all(x@parameters$fixed)){
    x@computed <- TRUE
    # Add fit:
    if (addfit){
      x <- addfit(x)
    }
    # FIXME: fis this Add MIs:
    # if (addMIs){
    #   x <- addMIs(x,analyticFisher=analyticFisher) 
    # }
    return(x)
  }
  
  optimizer <- x@optimizer
  # Default:
  if (optimizer == "default"){
    # if (x@model %in% c("varcov","lvm") && any(grepl("omega",x@parameters$matrix))){
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
    
    x@baseline_saturated$baseline <- runmodel(x@baseline_saturated$baseline, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE, addInformation = FALSE, analyticFisher = FALSE)
  }
  
  # Evaluate saturated model:
  
  if (!is.null(x@baseline_saturated$saturated) && !x@baseline_saturated$saturated@computed){
    if (verbose) message("Estimating saturated model...")
    # Run:
    x@baseline_saturated$saturated <- runmodel(x@baseline_saturated$saturated, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE, addInformation = FALSE, analyticFisher = FALSE)
  }
  
  
  
  # # nlminb control pars:
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
  #   optim.control <- modifyList(optim.control, control.nlminb)
  # }
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
  
  
  
  optim.control$par <- start
  optim.control$fn <- psychonetrics_fitfunction
  optim.control$model <- x
  if (level != "fitfunction"){
    optim.control$gr <- psychonetrics_gradient
  }
  
  # Add method:
  optim.control$method <- optimizer
  
  # Add bounds:
  if (optimizer %in% c("nlminb","L-BFGS-B","lbfgs")){

    optim.control$lower <- lower
    optim.control$upper <- upper
  }
  # Run model:
  curtry <- 1
  
  # If nlminb, add lavaan controls:
  if (optimizer == "nlminb"){
    if (is.null(optim.control$control)){
      optim.control$control<- list(eval.max=20000L,
                                     iter.max=10000L,
                                     trace=0L,
                                     #abs.tol=1e-20, ### important!! fx never negative
                                     abs.tol=(.Machine$double.eps * 10),
                                     # rel.tol=1e-10,
                                     rel.tol=1e-5,
                                     #step.min=2.2e-14, # in =< 0.5-12
                                     step.min=1.0, # 1.0 in < 0.5-21
                                     step.max=1.0,
                                     x.tol=1.5e-8,
                                     xf.tol=2.2e-14)
    }
  }
 
  
    repeat{
    tryres <- try({
      optim.out <- do.call(optimr,optim.control)
    }, silent = TRUE)    
    
    if (!is(tryres,"try-error") && !any(is.na(optim.out$par))){
      break
    } else {
      curtry <- curtry + 1
      if (curtry > maxtry){
        if (verbose){
          message("Model estimation failed and 'maxtry' reached. Returning error.")
          print(tryres)
        }        
        return(tryres)
      } else {
        if (verbose){
          message("Model estimation failed. Perturbing start values.")
        }        
        optim.control$par <- optim.control$par  + runif(length(optim.control$par),0, 0.1)
        optim.control$par[parMat(x) == "beta" | (rowMat(x) != colMat(x))] <- optim.control$par[parMat(x) == "beta" | (rowMat(x) != colMat(x))]/2
      }
    }
  }

  # optim.out <- do.call(optimr,optim.control)
  
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
  if (addInformation){
    if (verbose){
      message("Computing Fisher information...")
    }
    x@information <- psychonetrics_FisherInformation(x, analyticFisher)
    
    # if (verbose){
    #   message("Transpose...")
    # }
    if (!all(x@information == t(x@information))){
      x@information <- 0.5 * (x@information + t(x@information))
    }
    
    # if (verbose){
    #   message("Eigenvalues...")
    # }
    if (any(Re(eigen(x@information)$values) < -sqrt(.Machine$double.eps))){
      warning("Information matrix is not positive semi-definite. Model might not be identified.")
    }    
  }

  # }
  # Add fit:
  if (addfit){
    x <- addfit(x)
  }
  # Add MIs:
  if (addMIs){
    x <- addMIs(x,analyticFisher=analyticFisher, verbose = verbose) 
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