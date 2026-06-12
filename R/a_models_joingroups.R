allTheSame <- function(x)length(unique(x))==1

# Function to join models from different groups:
joingroups <- function(..., verbose = TRUE, log = TRUE, runmodel = FALSE,baseline_saturated = TRUE){
  # Store the groups:
  dots <- list(...)
  
  # Some checks:
  if (!all(sapply(dots,is,"psychonetrics"))){
    stop("All input to 'joingroups' must be of class 'psychonetrics;")
  }
  
  # Check if not single input:
  if (length(dots)==1){
    if (verbose) message("Only one model used as input to 'joingroups'. Nothing to do...")
    return(dots[[1]])
  }
  
  # If all models are not from the same class, stop:
  if (!allTheSame(sapply(dots,function(x)x@model))){
    stop ("Different models can not be combined.")
  }
  
  # Check if variables are the same:
  if (!allTheSame(sapply(dots,function(x)nrow(x@sample@variables)))){
    stop ("Models with different variables can not be combined.")
  }
  
  # Check more with the variables:
  vars <- sapply(dots,function(x)x@sample@variables$label)
  # sapply collapses to a vector for 1-variable models; keep matrix form
  # (variables in rows, models in columns):
  if (!is.matrix(vars)){
    vars <- matrix(vars, nrow = 1)
  }
  for (i in 1:nrow(vars)){
    if (!allTheSame(vars[i,])){
      stop ("Models with different variables can not be combined.")
    }
  }
  
  # Start with model 1 as base model:
  mod <- dots[[1]]
  mod@computed <- FALSE
  mod@equal <- "none"
  mod@baseline_saturated <- list()
  
  
  
  # combine the groups:
  mod@sample@groups <- do.call(rbind,lapply(dots,function(x)x@sample@groups))
  
  # Change id:
  mod@sample@groups$id <- seq_len(nrow(mod@sample@groups))
  
  # change the names (rename default single-group labels):
  defaultLabel <- mod@sample@groups$label %in% c("singlegroup","fullsample")
  mod@sample@groups$label[defaultLabel] <- paste0("group ",seq_len(sum(defaultLabel)))

  # Group labels must be unique (duplicate labels collapse group lookups by
  # name, e.g. in the C++ model matrix formation):
  mod@sample@groups$label <- make.unique(mod@sample@groups$label)
  
  # Change the sample statistics:
  mod@sample@covs <- do.call("c",lapply(dots,function(x)x@sample@covs))
  names(mod@sample@covs) <- mod@sample@groups$label
  
  mod@sample@cors <- do.call("c",lapply(dots,function(x)x@sample@cors))
  names(mod@sample@cors) <- mod@sample@groups$label
  
  mod@sample@means <- do.call("c",lapply(dots,function(x)x@sample@means))
  names(mod@sample@means) <- mod@sample@groups$label

  # Other per-group sample statistics (empty lists stay empty):
  mod@sample@thresholds <- do.call("c",lapply(dots,function(x)x@sample@thresholds))
  mod@sample@squares <- do.call("c",lapply(dots,function(x)x@sample@squares))
  mod@sample@fimldata <- do.call("c",lapply(dots,function(x)x@sample@fimldata))
  mod@sample@WLS.W <- do.call("c",lapply(dots,function(x)x@sample@WLS.W))
  mod@sample@WLS.Gamma <- do.call("c",lapply(dots,function(x)x@sample@WLS.Gamma))

  mod@sample@nobs <- sum(sapply(dots,function(x)x@sample@nobs))


  # Now add other models, add the parameters and the logs
  for (i in 2:length(dots)){
    newmod <- dots[[i]]@parameters
    newmod$group_id <- max(mod@parameters$group_id) + newmod$group_id
    newmod$par[newmod$par!=0] <- max(mod@parameters$par) + newmod$par[newmod$par!=0]

    mod@parameters <- rbind(mod@parameters,newmod)
    
    # obtain both logs:
    log1 <- mod@log
    log2 <- dots[[i]]@log
    
    # combined:
    comb <- c(log1,log2)
    
    # Dates:
    times <- sapply(comb,function(x)x@time)

    # Order the combined logs by time (c() and [ drop the S3 class, so
    # restore it before assigning to the slot):
    comb <- comb[order(times)]
    class(comb) <- "psychonetrics_log"
    mod@log <- comb
  }
  
  # Group name:
  # Group:
  mod@parameters$group <- mod@sample@groups$label[match(mod@parameters$group_id,mod@sample@groups$id)]
  
  # New model matrices:
  mod@modelmatrices <- formModelMatrices(mod)

  # Rebuild the parameter map matrix (M) for the combined parameter table;
  # leaving the old (single-group) M in place crashes the fit functions:
  if (mod@cpp){
    mod@extramatrices$M <- Mmatrix_cpp(mod@parameters)
  } else {
    mod@extramatrices$M <- Mmatrix(mod@parameters)
  }
  
  # Clear all estimates:
  mod@parameters$se[] <- NA
  mod@parameters$p[] <- NA
  mod@parameters$mi[] <- NA
  mod@parameters$pmi[] <- NA
  mod@parameters$mi_equal[] <- NA
  mod@parameters$pmi_equal[] <- NA
  
  # Add the baseline and saturated models:
  if (baseline_saturated){
    mod@baseline_saturated$baseline <- do.call(joingroups,c(lapply(dots,function(x)x@baseline_saturated$baseline),list(baseline_saturated=FALSE)))
    mod@baseline_saturated$saturated <- do.call(joingroups,c(lapply(dots,function(x)x@baseline_saturated$saturated),list(baseline_saturated=FALSE)))
  }

  # Write to log:
  if (log){
    # Add log:
    mod <- addLog(mod, paste0("Combined ",length(dots)," psychonetrics objects as groups!")) 
  }
  
  # Rerun if needed:
  if (runmodel){
    # NOTE: do NOT forward ... here; it contains the input models, which would
    # be positionally mismatched to runmodel's arguments.
    mod <- mod %>% runmodel(verbose=verbose)
  }
  
  mod
}