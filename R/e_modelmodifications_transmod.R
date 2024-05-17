transmod <- function(x,...,verbose,keep_computed = FALSE, log=TRUE, identify = TRUE){
  
  # Check if psychonetrics object:
  stopifnot(is(x, "psychonetrics"))
  
  # Set verbose:
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  # check if model is not scaled using variances:
  if (x@identification == "variance" && isTRUE(keep_computed)){
    stop("'keep_computed' cannot be used with identification = 'variance'; the transformed model needs to be run again.")
  }
  
  # Type of model:
  model <- x@model
  types <- x@types
  
  # dots:
  dots <- list(...)
  
  # length 1 and no name?
  if (is.null(names(dots)) && length(types) == length(dots)){
    names(dots) <- names(types)
  }
  
  # Else error if no names:
  if (is.null(names(dots)) ){
    stop("Named arguments representing the types to transform are required.")
  }
  
  # if dlvm1, fix two names:
  if (model == "dlvm1"){
    names(dots)[names(dots)=="between"] <- "between_latent"
    names(dots)[names(dots)=="within"] <- "within_latent"
  } else if (model == "varcov"){
    names(dots)[names(dots)=="type"] <- "y"
  }
  
  # Check if all names in types:
  if (any(!names(dots)%in%names(types))){
    stop(paste0("Not all types are valid types of the model ",model))
  }
  
  # Check if all types are not equal:
  eq <- TRUE
  for (t in seq_along(dots)){
    if (dots[[t]] != types[[names(dots)[t]]]){
      eq <- FALSE
    }
  }
  if (eq){
    message("No type to transform, returning the same model...")
    return(x)
  }
  
  # Else start transforming!!!
  
  # Type data frame:
  type_df <- rbind(
    
    # varcov:
    data.frame(
      family = "varcov",
      type = "y",
      appendix = ""
    ),
    
    # lvm:
    data.frame(
      family = "lvm",
      type = c("latent","residual"),
      appendix = c("_zeta","_epsilon")
    ),
    
    # dlvm1:
    data.frame(
      family = "dlvm1",
      type = c("within_latent","within_residual","between_latent","between_residual"),
      appendix = c("_zeta_within","_epsilon_within","_zeta_between","_epsilon_between")
    ),
    
    # ml_lvm:
    data.frame(
      family = "ml_lvm",
      type = c("within_latent","within_residual","between_latent","between_residual"),
      appendix = c("_zeta_within","_epsilon_within","_zeta_between","_epsilon_between")
    ),
    
    # var1:
    data.frame(
      family = "var1",
      type = c("zeta"),
      appendix = c("_zeta")
    ),
    
    # tsdlvm1:
    data.frame(
      family = "tsdlvm1",
      type = c("zeta","epsilon"),
      appendix = c("_zeta","_epsilon")
    ),
    
    # meta_varcov
    data.frame(
      family = "meta_varcov",
      type = c("y","randomEffects"),
      appendix = c("_y","_randomEffects")
    )
  )
  
  # Model matrices df:
  mod_matrices_df <- data.frame(
    type = c("cov","prec","chol","ggm","ggm","cor","cor"),
    matrix = c("sigma","kappa","lowertri","omega","delta","rho","SD")
  )
  
  # form the new model:
  newmod <- x
  
  # Update the matrices:
  newmod@modelmatrices <- formModelMatrices_cpp(newmod)
  
  # For each group:
  n_groups <- nrow(newmod@sample@groups)
  
  
  # For each type:
  for (t in seq_along(dots)){
    # appendix:
    appendix <- type_df$appendix[type_df$family == model & type_df$type == names(dots)[[t]]]
    
    # Obtain old and new type:
    old_type <- types[[names(dots)[t]]]
    new_type <- dots[[t]]
    
    # Obtain all current matrices #FIXME: I am lazy
    all_mats <- impliedcovstructures_cpp(newmod@modelmatrices,gsub("^\\_","",appendix),type=old_type,all=TRUE)
    
    # Obtain the relevant parameters:
    # Which parameters are these?
    if (appendix != ""){
      teststring <- appendix
    } else {
      teststring <- switch(x@types[[names(dots)[t]]],cov = "sigma", cor = "(rho)|(SD)", ggm = "(omega)|(delta)", prec = "kappa", chol = "lowertri"
      )
    }
    whichPars <- grepl(teststring,newmod@parameters$matrix)
    
    # FIXME: Not yet supported for models with equality constrains:
    if (any(duplicated(x@parameters$par[whichPars][x@parameters$par[whichPars]!=0]))){
      stop("'transmod' is not yet supported for models with equality constrains")
    }
    
    # Is the model diagonal or saturated?
    # Is the structure saturated?
    saturated <- all(newmod@parameters$est[whichPars]!=0)
    
    # Is the structure diagonal?
    diagonal <- all(newmod@parameters$est[whichPars & newmod@parameters$row == newmod@parameters$col]!=0) &
      all(newmod@parameters$est[whichPars & newmod@parameters$row != newmod@parameters$col]==0)
    
    # If not saturated or diagonal, stop:
    if (!(saturated || diagonal)){
      stop("Only diagonal and saturated model matrices can be transformed.")
    }
    
    
    
    # Overwrite the type:
    newmod@types[[names(dots)[t]]] <- new_type
    
    # Obtain covariance start values:
    exp_cov <- lapply(lapply(all_mats,"[[",paste0("sigma",appendix)), spectralshift)
    
    # variable names:
    varNames <- newmod@parameters$var1[whichPars & (newmod@parameters$row == newmod@parameters$col) & newmod@parameters$group_id == 1]
    
    # Obtain the new par table:
    if (saturated){

      new_par_tab <- do.call(generateAllParameterTables,matrixsetup_flexcov(sigma = "full",
                                                                            lowertri = "full",
                                                                            omega = "full",
                                                                            delta = "diag",
                                                                            kappa = "full",
                                                                            type = new_type,
                                                                            name= gsub("^\\_","",appendix),
                                                                            sampleStats= newmod@sample,
                                                                            equal = any(duplicated(newmod@parameters$par[whichPars & newmod@parameters$par > 0 & !newmod@parameters$fixed])),
                                                                            nNode = length(varNames),
                                                                            expCov = exp_cov,
                                                                            nGroup = nrow(newmod@sample@groups),
                                                                            labels = varNames,
                                                                            lassofix = FALSE))
    } else {
      new_par_tab <- do.call(generateAllParameterTables,matrixsetup_flexcov(sigma = "diag",
                                                                            lowertri = "diag",
                                                                            omega = "zero",
                                                                            delta = "diag",
                                                                            kappa = "diag",
                                                                            type = new_type,
                                                                            name= gsub("^\\_","",appendix),
                                                                            sampleStats= newmod@sample,
                                                                            equal = any(duplicated(newmod@parameters$par[whichPars & newmod@parameters$par > 0 & !newmod@parameters$fixed])),
                                                                            nNode = length(varNames),
                                                                            expCov = exp_cov,
                                                                            nGroup = nrow(newmod@sample@groups),
                                                                            labels = varNames,
                                                                            lassofix = FALSE))
    }
    
    
    # Remove underscore:
    if (appendix == ""){
      new_par_tab$partable$matrix <- gsub("\\_","",new_par_tab$partable$matrix)
      new_par_tab$mattable$name <- gsub("\\_","",new_par_tab$mattable$name)
    }
    
    # overwrite the parameter numbers:
    new_par_tab$partable$par <- max(x@parameters$par) + new_par_tab$partable$par
    
    # overwrite the parameters:
    newmod@parameters[whichPars,] <- new_par_tab$partable
    
    # Fix labels:
    newmod@parameters <- parRelabel(newmod@parameters)
    
    # overwrite the matrix:
    newmod@matrices <- rbind(
      newmod@matrices[!grepl(teststring,newmod@matrices$name),],
      new_par_tab$mattable
    )
    
    # Add log entry:
    if (log){
      newmod <- addLog(newmod,paste0("Transformed type ",names(dots)[t]," from ",old_type," to ",new_type))  
    }
    
    
  }
  
  # Add all matrices to the model:
  newmod@modelmatrices <- formModelMatrices_cpp(newmod)

  # Identify:
  if (identify){
    newmod <- identify(newmod)
  }
  
  # If the model was run, re-compute SEs:
  if (x@computed && x@identification != "variance"){
  
    # Add SEs:
    if (!all(is.na(x@parameters$se))){
      if (x@cpp){
        if (verbose){
          message("Adding standard errors...")
        }
        newmod <- addSEs_cpp(newmod,verbose=verbose)  
      } else {
        newmod <- addSEs(newmod, verbose=verbose)
      }
    }
    
    
    
    
    # Add MIs:
    if (!all(is.na(x@parameters$mi))){
      newmod <- addMIs(newmod)
    }
    
    
  }
  
  # Keep computed?
  if (!keep_computed){
    newmod@computed <- FALSE
  }
  

  
  # Return the new model:
  return(newmod)
}
