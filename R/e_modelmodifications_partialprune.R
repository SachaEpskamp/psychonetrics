partialprune <- function(
  x, # Model
  alpha = 0.01, # Significance
  matrices, # Automatically chosen
  verbose,
  combinefun = unionmodel,
  return = c("best","partialprune","union_equal","prune"),
  criterion = "bic",
  best = c("lowest","highest"),
  ...){
  # Verbose:
  if (missing(verbose)){
    verbose <- x@verbose
  }
  fixed <- NULL
  mi_free <- NULL
  
  # If not run, run model:
  if (!x@computed){
    x <- x %>% runmodel(..., verbose = verbose)
  }
  
  # Test criterion:
  if (!criterion %in% names( x@fitmeasures)){
    stop("'criterion' is not a valid psychonetrics fit measure \ information criterion")
  }
  
  # Return argument:
  return <- match.arg(return)
  best <- match.arg(best)
  
  # Matrices:
  if (missing(matrices)){
    if (x@model == "varcov"){
      if (x@submodel == "ggm"){
        matrices <- "omega"
      } else if (x@submodel == "prec"){
        matrices <- "kappa"
      } else if (x@submodel == "chol"){
        matrices <- "lowertri"
      } else if (x@submodel == "cov"){
        matrices <- "sigma"
      } 
    }  else if (x@model == "meta_varcov"){
      
      matrices <- character(0)
      if (x@types$y == "ggm"){
        matrices <- c(matrices,"omega_y")
      } else if (x@types$y == "prec"){
        matrices <- c(matrices,"kappa_y") 
      } 
      
      # Random effects:
      if (x@types$randomEffects == "ggm"){
        matrices <- c(matrices,"omega_randomEffects")
      } else if (x@types$randomEffects == "prec"){
        matrices <- c(matrices,"kappa_randomEffects") 
      } 
      
      
    }  else if (x@model == "lvm"){
      
      # Only add GGM structures to search:
      matrices <- character(0)
      if (x@types$latent == "ggm"){
        matrices <- c(matrices,"omega_zeta")
      } else if (x@types$latent == "prec"){
        matrices <- c(matrices,"kappa_zeta")
      }
      
      if (x@types$residual == "ggm"){
        matrices <- c(matrices,"omega_epsilon")
      } else if (x@types$residual == "prec"){
        matrices <- c(matrices,"kappa_epsilon")
      }
      
      
    } else if (x@model == "lnm"){
      matrices <- "omega_eta"
    } else if (x@model == "rnm"){
      matrices <- "omega_epsilon"
    }  else if (x@model == "gvar"){
      matrices <- c("beta","omega_zeta")
    }   else if (x@model == "var1"){
      matrices <- c("beta")
      if (x@types$zeta == "prec"){
        matrices <- c(matrices,"kappa_zeta")
      } else if (x@types$zeta == "ggm"){
        matrices <- c(matrices,"omega_zeta")
      }
      
    } else if (x@model == "panelvar1"){
      matrices <- c("beta")
      if (x@types$contemporaneous == "prec"){
        matrices <- c(matrices,"kappa_zeta")
      } else if (x@types$contemporaneous == "ggm"){
        matrices <- c(matrices,"omega_zeta")
      }
      
      if (x@types$between == "prec"){
        matrices <- c(matrices,"kappa_mu")
      } else if (x@types$between == "ggm"){
        matrices <- c(matrices,"omega_mu")
      }
      
    } else if (x@model %in% c("ml_lvm","dlvm1")){
      if (x@model == "dlvm1") {
        matrices <- c("beta")
      } else {
        matrices <- character(0)
      }
      if (x@types$within_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_within")
      } else if (x@types$within_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_within")
      }
      
      if (x@types$within_residual == "prec"){
        matrices <- c(matrices,"kappa_epsilon_within")
      } else if (x@types$within_residual == "ggm"){
        matrices <- c(matrices,"omega_epsilon_within")
      }
      
      if (x@types$between_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_between")
      } else if (x@types$between_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_between")
      }
      
      if (x@types$between_residual == "prec"){
        matrices <- c(matrices,"kappa_epsilon_between")
      } else if (x@types$between_residual == "ggm"){
        matrices <- c(matrices,"omega_epsilon_between")
      }
      
    } else if (x@model == "tsdlvm1"){
      matrices <- c("beta")
      if (x@types$zeta == "prec"){
        matrices <- c(matrices,"kappa_zeta")
      } else if (x@types$zeta == "ggm"){
        matrices <- c(matrices,"omega_zeta")
      }
      
      if (x@types$epsilon == "prec"){
        matrices <- c(matrices,"kappa_epsilon")
      } else if (x@types$epsilon == "ggm"){
        matrices <- c(matrices,"omega_epsilon")
      }
      
    }  else if (x@model == "Ising"){
      matrices <- c("omega")
      
    }  else stop("No default argument for 'matrices' for current model.")
  }
  
  # Prune first:
  if (verbose) message("Pruning model...")
  mod_prune <- prune(x,alpha=alpha,verbose=FALSE,runmodel=TRUE,matrices=matrices,...) 
  
  # if not empty, look for equality:
  if (!all(mod_prune@parameters$est[mod_prune@parameters$matrix%in%matrices&!mod_prune@parameters$fixed]==0)){
    
    
    # Combine models:
    if (verbose) message("Combining models...")
    
    # First union or intersection:
    mod_union <- combinefun(mod_prune, matrices = matrices)
    
    # Then set equal:
    for (m in seq_along(matrices)){
      mod_union <- groupequal(mod_union,matrices[m], verbose = FALSE)  
    }
  
    # Then run:
    mod_union <- runmodel(mod_union, verbose = FALSE)
    
    if (verbose) message("Partial pruning...")
    curMod <- mod_union
    
    
    repeat{
      pars <- curMod@parameters
      
      # Make a data frame in which the equality free parameters are summed:
      miDF <- pars %>% filter(!fixed) %>% group_by(.data[["row"]],.data[["col"]],.data[["matrix"]]) %>%
        filter(matrix %in% matrices) %>%
        summarize(mi_free = sum(mi_free)) %>% 
        arrange(-mi_free)
      
      # break if empty:
      if (nrow(miDF)==0){
        break
      }
      
      # Free the best parameter:
      propMod <- curMod %>% 
        groupfree(miDF$matrix[1],miDF$row[1],miDF$col[1]) %>% 
        runmodel
      
      # Test BIC:
      if (propMod@fitmeasures$bic < curMod@fitmeasures$bic){
        curMod <- propMod
      } else {
        break
      }
    }
    mod_partialpooled <- curMod %>% prune(alpha=alpha,verbose=FALSE,runmodel=TRUE,matrices=matrices,...) 
    
    # Select best model:
    mods <- list(x, mod_prune, mod_union, mod_partialpooled)
    
    if (best == "lowest"){
      best <- which.min(c(
        x@fitmeasures[[criterion]],
        mod_prune@fitmeasures[[criterion]],
        mod_union@fitmeasures[[criterion]],
        mod_partialpooled@fitmeasures[[criterion]]
      ))
    } else {
      best <- which.max(c(
        x@fitmeasures[[criterion]],
        mod_prune@fitmeasures[[criterion]],
        mod_union@fitmeasures[[criterion]],
        mod_partialpooled@fitmeasures[[criterion]]
      ))
      
    }
    
    bestmod <- mods[[best]]
    
    # Which model to return?
    if (return == "best"){
      return(bestmod)  
    } else if (return == "partialprune") {
      return(mod_partialpooled)  
    } else if (return == "union_equal") {
      return(mod_union)  
    } else if (return == "prune") {
      return(mod_prune)  
    } else {
      stop("Incorrect 'return' argument.")
    }
    
    
  } else {
    
    if (verbose) message("Model matrices are empty after pruning, returning empty matrices")
    return(mod_prune)
  }

  

}