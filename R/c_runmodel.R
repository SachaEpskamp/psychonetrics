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
    optim.control,
    # maxtry = 5,
    analyticFisher = TRUE,
    warn_improper = FALSE,
    warn_gradient = TRUE,
    warn_bounds = TRUE,
    return_improper = TRUE,
    bounded = TRUE,
    approximate_SEs=FALSE
    # cholesky_start # If TRUE, a model is formed with Cholesky decompositions first which is run for obtaining starting values.
){
  # cholesky_start <- FALSE
  # Cholesky start:
  # if (missing(cholesky_start)){
  #   # # Don't do this if the model was evaluated earlier:
  #   # cholesky_start <- x@distribution == "Gaussian" && !any(sapply(x@log,slot,"event")=="Evaluated model")
  # 
  #   # Only Gaussian:
  #   cholesky_start <- x@distribution == "Gaussian"
  # 
  # }
  # if (!is.logical(cholesky_start)){
  #   stop("'cholesky_start' must be a logical argument")
  # }
  # if (isTRUE(cholesky_start) && x@distribution != "Gaussian"){
  #   stop("'cholesky_start' is only supported for Gaussian models")
  # }
  
  if (!missing(optim.control)){
    warning("'optim.control' is deprecated and will be removed in a future version. Please use setoptimizer(..., optim.args = ...).")
    x@optim.args <- optim.control
  }
  
  # Set optim args:
  optim.control <- x@optim.args
  
  
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
  # 
  
  
  ### CHOLESKY DECOMPOSITION STARTING VALUES ###
  # if (cholesky_start){
  #   
  #   # FIXME: a data frame with all types and classes in psychonetrics:
  #   type_df <- rbind(
  #     
  #     # varcov:
  #     data.frame(
  #       family = "varcov",
  #       type = "y",
  #       appendix = ""
  #     ),
  #     
  #     # lvm:
  #     data.frame(
  #       family = "lvm",
  #       type = c("latent","residual"),
  #       appendix = c("_zeta","_epsilon")
  #     ),
  #     
  #     # dlvm1:
  #     data.frame(
  #       family = "dlvm1",
  #       type = c("within_latent","within_residual","between_latent","between_residual"),
  #       appendix = c("_zeta_within","_epsilon_within","_zeta_between","_epsilon_between")
  #     ),
  #     
  #     # ml_lvm:
  #     data.frame(
  #       family = "ml_lvm",
  #       type = c("within_latent","within_residual","between_latent","between_residual"),
  #       appendix = c("_zeta_within","_epsilon_within","_zeta_between","_epsilon_between")
  #     ),
  #     
  #     # var1:
  #     data.frame(
  #       family = "var1",
  #       type = c("zeta"),
  #       appendix = c("_zeta")
  #     ),
  #     
  #     # tsdlvm1:
  #     data.frame(
  #       family = "tsdlvm1",
  #       type = c("zeta","epsilon"),
  #       appendix = c("_zeta","_epsilon")
  #     ),
  #     
  #     # meta_varcov
  #     data.frame(
  #       family = "meta_varcov",
  #       type = c("y","randomEffects"),
  #       appendix = c("_y","_randomEffects")
  #     )
  #   )
  #   
  #   
  #   
  #   # Form the Cholesky model:
  #   chol_mod <- x
  #   
  #   # Update the matrices:
  #   chol_mod@modelmatrices <- formModelMatrices_cpp(chol_mod)
  #   
  #   # Check if any type is not chol?
  #   curtypes <- unlist(x@types)
  #   
  #   # New mod matrices to be stored:
  #   new_modmatrices <- list()
  #   
  #   # look over all types that are not chol (if all are chol, do nothing):
  #   for (type in which(curtypes!="chol")){
  #     
  #     # What is the appendix?
  #     appendix <- type_df$appendix[type_df$family == x@model & type_df$type == names(curtypes!="chol")[type]]
  #     
  #     # Test string:
  #     
  #     # Which parameters are these?
  #     if (appendix != ""){
  #       teststring <- appendix
  #     } else {
  #       teststring <- switch(x@types[[type]],cov = "sigma", cor = "(rho)|(SD)", ggm = "(omega)|(delta)", prec = "kappa", chol = "lowertri"
  #       )
  #     }
  #     whichPars <- grepl(teststring,chol_mod@parameters$matrix)
  # 
  #     # Is the structure saturated?
  #     saturated <- all(!chol_mod@parameters$fixed[whichPars])
  #     
  #     # Is the structure diagonal?
  #     diagonal <- all(!chol_mod@parameters$fixed[whichPars & chol_mod@parameters$row == chol_mod@parameters$col]) &
  #       all(chol_mod@parameters$fixed[whichPars & chol_mod@parameters$row != chol_mod@parameters$col]) 
  #     
  #     # FIXME: Only do this if the model is saturated or diagonal:
  #     if (saturated || diagonal){
  #       
  #       # Make the type a cholesky decomposition:
  #       chol_mod@types[[type]] <- "chol"
  #       
  #       # Extract the covariance matrix:
  #       # FIXME: Only sigma is needed, but I am lazy
  #       all_mats <- impliedcovstructures_cpp(chol_mod@modelmatrices,gsub("^\\_","",appendix),type=x@types[[type]],all=FALSE)
  #       
  #       # Obtain covariance start values:
  #       exp_cov <- lapply(lapply(all_mats,"[[",paste0("sigma",appendix)), spectralshift)
  #       
  #       # variable names:
  #       varNames <- chol_mod@parameters$var1[whichPars & (chol_mod@parameters$row == chol_mod@parameters$col) & chol_mod@parameters$group_id == 1]
  # 
  #       # Obtain the new par table:
  #       if (saturated){
  #         new_par_tab <- generateAllParameterTables(matrixsetup_flexcov(sigma = "full",
  #                                                                       lowertri = "full",
  #                                                                       omega = "full",
  #                                                                       delta = "diag",
  #                                                                       kappa = "full",
  #                                                                       type = "chol",
  #                                                                       name= gsub("^\\_","",appendix),
  #                                                                       sampleStats= chol_mod@sample,
  #                                                                       equal = any(duplicated(chol_mod@parameters$par[whichPars & chol_mod@parameters$par > 0 & !chol_mod@parameters$fixed])),
  #                                                                       nNode = length(varNames),
  #                                                                       expCov = exp_cov,
  #                                                                       nGroup = nrow(chol_mod@sample@groups),
  #                                                                       labels = varNames)[[1]])
  #       } else {
  #         new_par_tab <- generateAllParameterTables(matrixsetup_flexcov(sigma = "diag",
  #                                                                       lowertri = "diag",
  #                                                                       omega = "zero",
  #                                                                       delta = "diag",
  #                                                                       kappa = "diag",
  #                                                                       type = "chol",
  #                                                                       name= gsub("^\\_","",appendix),
  #                                                                       sampleStats= chol_mod@sample,
  #                                                                       equal = any(duplicated(chol_mod@parameters$par[whichPars & chol_mod@parameters$par > 0 & !chol_mod@parameters$fixed])),
  #                                                                       nNode = length(varNames),
  #                                                                       expCov = exp_cov,
  #                                                                       nGroup = nrow(chol_mod@sample@groups),
  #                                                                       labels = varNames)[[1]])
  #       }
  #       
  #       # Remove underscore:
  #       if (appendix == ""){
  #         new_par_tab$partable$matrix <- gsub("\\_","",new_par_tab$partable$matrix)
  #         new_par_tab$mattable$name <- gsub("\\_","",new_par_tab$mattable$name)
  #       }
  #       
  #       # overwrite the parameter numbers:
  #       new_par_tab$partable$par <- chol_mod@parameters$par[whichPars]
  #       
  #       # overwrite the parameters:
  #       chol_mod@parameters[whichPars,] <- new_par_tab$partable
  #       
  #       # overwrite the matrix:
  #       chol_mod@matrices <- rbind(
  #         chol_mod@matrices[!grepl(teststring,chol_mod@matrices$name),],
  #         new_par_tab$mattable
  #       )
  #       
  #     }
  #   }
  # 
  #   
  #   # Remove baseline and saturated:
  #   chol_mod@baseline_saturated$baseline <- NULL
  #   chol_mod@baseline_saturated$saturated <- NULL
  #   
  #   if (verbose){
  #     message("Running initial Cholesky decomposition model for starting values...")
  #   }
  #   
  #   # After including Cholesky decompositions, run the model:
  #   tryres_chol <- try({
  #     suppressMessages({
  #       suppressWarnings({
  #         chol_mod <- runmodel(chol_mod, addfit = FALSE, addMIs = FALSE, verbose = FALSE,addSEs=FALSE, addInformation = FALSE, analyticFisher = FALSE, cholesky_start = FALSE)
  #       })
  #     })
  #   })
  #   
  #   # Only if this succeeded overwrite the starting values..
  #   if (!is(tryres_chol,"try-error")){
  #     
  #     # first overwrite most start values, this will handle intercepts, loadings, etcetera:
  #     x@parameters$par[!x@parameters$fixed] <- chol_mod@parameters$par[!x@parameters$fixed]
  #     
  #     # FIXME: Code repetition below...
  #     
  #     # Then we need to loop again...
  #     for (type in which(curtypes!="chol")){
  #       
  #       # What is the appendix?
  #       appendix <- type_df$appendix[type_df$family == x@model & type_df$type == names(curtypes!="chol")[type]]
  #       
  #       # Which parameters are these?
  #       if (appendix != ""){
  #         teststring <- appendix
  #       } else {
  #         teststring <- switch(x@types[[type]],cov = "sigma", cor = "(rho)|(SD)", ggm = "(omega)|(delta)", prec = "kappa", chol = "lowertri"
  #         )
  #       }
  #       
  #       whichPars <- grepl(teststring,x@parameters$matrix)
  #       
  #       # Is the structure saturated?
  #       saturated <- all(!chol_mod@parameters$fixed[whichPars])
  #       
  #       # Is the structure diagonal?
  #       diagonal <- all(!chol_mod@parameters$fixed[whichPars & chol_mod@parameters$row == chol_mod@parameters$col]) &
  #         all(chol_mod@parameters$fixed[whichPars & chol_mod@parameters$row != chol_mod@parameters$col]) 
  #       
  #       # FIXME: Only do this if the model is saturated or diagonal:
  #       if (saturated || diagonal){
  #         
  #         # FIXME: Only the relevant matrices are needed, but I am lazy
  #         all_mats <- impliedcovstructures_cpp(chol_mod@modelmatrices,gsub("^\\_","",appendix),type="chol",all=FALSE)
  # 
  #         # Obtain covariance start values:
  #         exp_cov <- lapply(lapply(all_mats,"[[",paste0("sigma",appendix)), spectralshift)
  #         
  #         # variable names:
  #         varNames <- chol_mod@parameters$var1[whichPars & chol_mod@parameters$row == chol_mod@parameters$col & chol_mod@parameters$group_id == 1]
  #         
  #         # Obtain the new par table:
  #         new_par_tab <- do.call(generateAllParameterTables,matrixsetup_flexcov(sigma = "full",
  #                                                                       lowertri = "full",
  #                                                                       omega = "full",
  #                                                                       delta = "diag",
  #                                                                       kappa = "full",
  #                                                                       type = x@types[[type]],
  #                                                                       name= gsub("^\\_","",appendix),
  #                                                                       sampleStats= chol_mod@sample,
  #                                                                       equal = FALSE,
  #                                                                       nNode = length(varNames),
  #                                                                       expCov = exp_cov,
  #                                                                       nGroup = nrow(chol_mod@sample@groups),
  #                                                                       labels = varNames))
  #        
  #         # Overwrite the starting values:
  #         x@parameters$est[whichPars & !x@parameters$fixed] <- new_par_tab$partable$est[!x@parameters$fixed[whichPars]]
  #       
  #         # For the omega and rho matrices, bound parameters between -0.5 and 0.5:
  #         if (x@types[[type]] %in% c("ggm","cor")){
  #           edges <- whichPars & !x@parameters$fixed & grepl("(omega)|(rho)",x@parameters$matrix)
  #           x@parameters$est[edges] <- pmax(pmin(x@parameters$est[edges], 0.5),-0.5)
  #         }
  #         
  #       }
  #     }
  #   }
  # }

  
  
  # FIXME: Ugly loop to check for start values
  # trystart <- 1
  # while (trystart < 3){
  # Start and bounds:
  start <- parVector(x)
  lower <- lowerBound(x)
  upper <- upperBound(x)
  
  oldstart <- start
  
  if (verbose) message("Estimating model...")
  # Form optimizer arguments:
  
  
  if (grepl("cpp",optimizer)){
    
    suppressWarnings({
      tryres <- try({
        x <- psychonetrics_optimizer(x, lower, upper, gsub("cpp_","",optimizer), bounded)
      }, silent = TRUE) 
    })
    
    if (is(tryres,"try-error") && !any(is.na(parVector(x)))){
      
      suppressWarnings({
        tryres2 <- try({
          # browser()
          x <- updateModel(oldstart, x)
          x <- psychonetrics_optimizer(emergencystart(x), lower, upper, gsub("cpp_","",optimizer))
        }, silent = TRUE)    
      })
      
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
    if (optimizer %in% c("nlminb","L-BFGS-B","lbfgs") && bounded){
      
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
                                     abs.tol=sqrt(.Machine$double.eps),
                                     rel.tol=sqrt(.Machine$double.eps),
                                     step.min=1.0,
                                     step.max=1.0,
                                     x.tol=1.5e-8,
                                     xf.tol=2.2e-14)
        # 
        # optim.control$control<- list(maxfeval=20000L,
        #                              maxit=10000L,
        #                              trace=0L)
      }
    }
    
    suppressWarnings({
      tryres <- try({
        optim.out <- do.call(optimr_fake,optim.control)
      }, silent = TRUE)    
    })
    
    if (is(tryres,"try-error") || any(is.na(optim.out$par))){
      # Try with emergencystart:
      x <- updateModel(oldstart, x)
      optim.control$par <- parVector(emergencystart(x))
      
      suppressWarnings({
        tryres2 <- try({
          optim.out <- do.call(optimr_fake,optim.control)
        }, silent = TRUE)    
      })
      
      # If still an error, break:
      if (is(tryres2,"try-error") || any(is.na(optim.out$par))){
        stop(paste("Model estimation resulted in an error that could not be recovered. Try using a different optimizer with setoptimizer(...) or using different starting values. Error:\n\n",tryres2))
        
      }
      
    }
    
    
    optimresults <- optim.out
    optimresults$optimizer <- optimizer
    class(optimresults) <- c("list",class(optimresults))
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
    # 
    # 
    # # If information is not positive definite, try to fix:
    # if (!sympd_cpp(x@information, semi = FALSE)){
    #   
    #   warning("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try standardizing data or adjusting starting values.")    
    #   
    # }
    
    if (!proper && warn_improper){
      warning("The optimizer encountered at least one non-positive definite matrix and used an approximate inverse in parameter estimation. This is likely not a problem, but make sure to inspect if the parameters look reasonable before evaluating model fit.")
    }
    # 
    #   # Try to recover using start values:
    #   # if (trystart == 1){
    #   #   trystart <- 2
    #   #   x <- updateModel(oldstart, x)
    #   #   x <- emergencystart(x)
    #   # } else {
    #     trystart <- 3
    #     if (return_improper){
    #       if (warn_improper){
    #         warning("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try standardizing data or adjusting starting values.")    
    #       }
    #     } else {
    #       stop("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try standardizing data or adjusting starting values.")
    #     }
    #   # }
    #   
    #   
    # } else if (!proper){
    #   # Check for a non-proper computation:
    #   
    #   # Either return an error or a warning:
    #   if (!return_improper){
    #     stop("The optimizer encountered at least one non-positive definite matrix and used a pseudoinverse in parameter estimation. To return results anyway, set return_improper = TRUE.")
    #   } else{
    #     if (warn_improper){
    #       warning("The optimizer encountered at least one non-positive definite matrix and used a pseudoinverse in parameter estimation. Inspect if the parameters look reasonable before evaluating model fit.")
    #       trystart <- 3
    #     }
    #     
    #   }
    #   
    #   
    # } else {
    #   trystart <- 3
    # }
    # 
    # 
  }
  
  # Check bounds:
  if (bounded && warn_bounds){
    if (!all(x@parameters$fixed)){
      if (any(x@parameters$est[!x@parameters$fixed] <= x@parameters$minimum[!x@parameters$fixed] + sqrt(.Machine$double.eps)) ||
          any(x@parameters$est[!x@parameters$fixed] >= x@parameters$maximum[!x@parameters$fixed] - sqrt(.Machine$double.eps))){
        
        # which_wrong <- which(x@parameters$est[!x@parameters$fixed] <= x@parameters$lower[!x@parameters$fixed] + sqrt(.Machine$double.eps) |
        #                        x@parameters$est[!x@parameters$fixed] >= x@parameters$upper[!x@parameters$fixed] - sqrt(.Machine$double.eps))
        
        warning("One or more parameters were estimated to be near its bounds. This may be indicative of, for example, a Heywood case, but also of an optimization problem. Interpret results and fit with great care. For unconstrained estimation, set bounded = FALSE.")
        
      }
      
    }
    
  }
  # 
  # if (!sympd_cpp(x@information) || (!return_improper && !proper)){
  #   if (trystart == 1){
  #     trystart <- 2
  #     x <- updateModel(oldstart, x)
  #     x <- emergencystart(x)
  #   } else {
  #     trystart <- 3
  #     if (return_improper){
  #       if (warn_improper){
  #         warning("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).")    
  #       }
  #     } else {
  #       stop("Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).")
  #     }
  #   }
  #   
  # } else {
  #   trystart <- 3
  # }
  # if (any(Re(eigen(x@information)$values) < -sqrt(.Machine$double.eps))){
  #   warning("Information matrix is not positive semi-definite. Model might not be identified.")
  # }    
  
  #   } else {
  #     trystart <- 3
  #   }
  # }
  # 
  
  # }
  # Add baseline or saturated if needed:
  if (isBaseline){
    x@baseline_saturated$baseline <- x
  }
  if (isSaturated){
    x@baseline_saturated$saturated <- x
  }
  
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
      x <- addSEs_cpp(x,verbose=verbose,approximate_SEs=approximate_SEs)  
    } else {
      x <- addSEs(x, verbose=verbose,approximate_SEs=approximate_SEs)
    }
    
  }
  
  
  
  if (log){
    # Add log:
    x <- addLog(x, "Evaluated model")    
  }
  
  # Warn about the gradient?
  if (warn_gradient){
    
    if (x@cpp){
      grad <- psychonetrics_gradient_cpp(parVector(x),x)
    } else {
      grad <- psychonetrics_gradient(parVector(x),x)
    }
    
    
    if (mean(abs(grad)) > 1){
      warning("Model might not have converged properly: mean(abs(gradient)) > 1.")
    }
  }
  
  
  # Return model:
  return(x)
}