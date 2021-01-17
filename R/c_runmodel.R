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
  verbose,
  # optimizer = c("default","ucminf","nlminb"),
  optim.control = list(),
  # maxtry = 5,
  analyticFisher = TRUE,
  return_improper = FALSE
  # inverseHessian = TRUE
){
  if (missing(verbose)){
    verbose <- x@verbose
  }
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
    if (x@distribution == "Gaussian" && any(grepl("omega",x@parameters$matrix))){
      optimizer <- "nlminb"
    } else {
      optimizer <- "ucminf"
    }
    # optimizer <- "psychonetrics_BFGS"
  }
  
  
  # optimizer <- match.arg(optimizer)
  level <- match.arg(level)
  if (!is(x,"psychonetrics")){
    stop("input is not a 'psychonetrics' object")
  }
  
  
  if (!is.null(x@baseline_saturated$baseline)){
    # Check if model happens to be baseline model:
    isBaseline <- identical(x@baseline_saturated$baseline@parameters$par, x@parameters$par)    
  } else {
    isBaseline <- FALSE
  }
  
  if (!is.null(x@baseline_saturated$saturated)){
    # Check if model happens to be saturated model:
    isSaturated <- identical(x@baseline_saturated$saturated@parameters$par, x@parameters$par)
  } else {
    isSaturated <- FALSE
  }
  
  # Evaluate baseline model:
  if (!isBaseline && !is.null(x@baseline_saturated$baseline) && !x@baseline_saturated$baseline@computed){
    if (verbose) message("Estimating baseline model...")
    # Run:
    x@baseline_saturated$baseline@optimizer <- optimizer
    x@baseline_saturated$baseline <- runmodel(x@baseline_saturated$baseline, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE, addInformation = FALSE, analyticFisher = FALSE)
  }
  
  # Evaluate saturated model:
  
  if (!isSaturated && !is.null(x@baseline_saturated$saturated) && !x@baseline_saturated$saturated@computed){
    if (verbose) message("Estimating saturated model...")
    x@baseline_saturated$saturated@optimizer <- optimizer
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
  
  # FIXME: Ugly loop to check for start values
  trystart <- 1
  while (trystart < 3){
    # Start and bounds:
    start <- parVector(x)
    lower <- lowerBound(x)
    upper <- upperBound(x)
    
    oldstart <- start

    if (verbose) message("Estimating model...")
    # Form optimizer arguments:
    
    
    if (grepl("cpp",optimizer)){
      
      
      
      
      
      tryres <- try({
        x <- psychonetrics_optimizer(x, lower, upper, gsub("cpp_","",optimizer))
      }, silent = TRUE)    
      
      if (is(tryres,"try-error") && !any(is.na(parVector(x)))){
        
        tryres2 <- try({
          # browser()
          x <- updateModel(oldstart, x)
          x <- psychonetrics_optimizer(emergencystart(x), lower, upper, gsub("cpp_","",optimizer))
        }, silent = TRUE)    
        
        # If still an error, break:
        if (is(tryres2,"try-error") && !any(is.na(parVector(x)))){
          stop(paste("Model estimation resulted in an error that could not be recovered. Try using a different optimizer with setoptimizer(...) or using different starting values. Error:\n\n",tryres2))
          
        }
        
      }
      
      
      x@objective <- x@optim$value
      
      x <- updateModel(parVector(x),x,updateMatrices = TRUE) # FIXME: Move this to C++!
      
      
      
      
      
      
      
    } else {
      optim.control$par <- start
      if (x@cpp){
        optim.control$fn <- psychonetrics_fitfunction_cpp
      } else {
        optim.control$fn <- psychonetrics_fitfunction
      }
      
      optim.control$model <- x
      if (level != "fitfunction"){
        if (x@cpp){
          optim.control$gr <- psychonetrics_gradient_cpp
        } else {
          optim.control$gr <- psychonetrics_gradient
        }
        
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
          # 
          # optim.control$control<- list(maxfeval=20000L,
          #                              maxit=10000L,
          #                              trace=0L)
        }
      }
      
      
      tryres <- try({
        optim.out <- do.call(optimr_fake,optim.control)
      }, silent = TRUE)    

      if (is(tryres,"try-error") || any(is.na(optim.out$par))){
        # Try with emergencystart:
        x <- updateModel(oldstart, x)
        optim.control$par <- parVector(emergencystart(x))
        
        tryres2 <- try({
          optim.out <- do.call(optimr_fake,optim.control)
        }, silent = TRUE)    
        
        # If still an error, break:
        if (is(tryres2,"try-error") || any(is.na(optim.out$par))){
          stop(paste("Model estimation resulted in an error that could not be recovered. Try using a different optimizer with setoptimizer(...) or using different starting values. Error:\n\n",tryres2))
          
        }
        
      }
      
    
      optimresults <- optim.out
      optimresults$optimizer <- optimizer
      x@optim <- optimresults
      # x@computed <- TRUE
      x@objective <- optimresults$value
      x <- updateModel(optim.out$par,x,updateMatrices = TRUE)
    }
    
    
    x@optimizer <- optimizer
    
    # optim.out <- do.call(optimr,optim.control)
    
    # Update model:
    
    
    # Make list:
    # optimresults <- list(
    #   par = optim.out$par,
    #   value = optim.out$value,
    #   message = optim.out$message,
    #   optimizer = optimizer
    # )
    
    
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
    
    # x@optim <- optimresults
    x@computed <- TRUE
    # x@objective <- optimresults$value
    
    # Add information:
    # if (!is.null(x@fitfunctions$information)){
    if (addInformation){
      if (verbose){
        message("Computing Fisher information...")
      }
      
      if (x@cpp){
        x@information <- psychonetrics_FisherInformation_cpp(x, analyticFisher)
      } else {
        x@information <- psychonetrics_FisherInformation(x, analyticFisher)
      }
      
      
      # if (verbose){
      #   message("Transpose...")
      # }
      # if (!all(x@information == t(x@information))){
      #   x@information <- 0.5 * (x@information + t(x@information))
      # }
      # 
      # if (verbose){
      #   message("Eigenvalues...")
      # }
      # Check if analysis was proper:
      proper <- TRUE
      if (!is.null(x@modelmatrices[[1]]$proper)){
        propers <- sapply(x@modelmatrices,"[[","proper")
        proper <- all(propers)
      }
      
      if (!sympd_cpp(x@information) || (!return_improper && !proper)){
        if (trystart == 1){
          trystart <- 2
          x <- updateModel(oldstart, x)
          x <- emergencystart(x)
        } else {
          trystart <- 3
          if (return_improper){
            warning("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).")  
          } else {
            stop("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).")
          }
        }
        
      } else {
        trystart <- 3
      }
      # if (any(Re(eigen(x@information)$values) < -sqrt(.Machine$double.eps))){
      #   warning("Information matrix is not positive semi-definite. Model might not be identified.")
      # }    
      
    } else {
      trystart <- 3
    }
  }
  
  
  # }
  # Add fit:
  if (addfit){
    x <- addfit(x, verbose=verbose)
  }
  # Add MIs:
  if (addMIs){
    x <- addMIs(x,analyticFisher=analyticFisher, verbose = verbose) 
  }
  # Add SEs:
  if (addSEs){
    if (x@cpp){
      if (verbose){
        message("Adding standard errors...")
      }
      x <- addSEs_cpp(x)  
    } else {
      x <- addSEs(x, verbose=verbose)
    }
    
  }
  
  # Add baseline or saturated if needed:
  if (isBaseline){
    x@baseline_saturated$baseline <- x
  }
  if (isSaturated){
    x@baseline_saturated$saturated <- x
  }
  
  if (log){
    # Add log:
    x <- addLog(x, "Evaluated model")    
  }

  
  
  # Return model:
  return(x)
}