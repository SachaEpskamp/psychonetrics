partialprune <- function(
  x, # Model
  alpha = 0.01, # Significance
  matrices, # Automatically chosen
  verbose,
  combinefun = c("unionmodel","intersectionmodel","identity"), # Allowed options: unionmodel, intersectionmodel, or identity
  return = c("partialprune","best","union_equal","prune"),
  criterion = "bic",
  best = c("lowest","highest"),
  final_prune = c("saturated","partialprune"),
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
  final_prune <- match.arg(final_prune)
  
  if (identical(combinefun,unionmodel)){
    combinefun <- "unionmodel"
  }
  if (identical(combinefun,intersectionmodel)){
    combinefun <- "intersectionmodel"
  }
  
  combinefun <- match.arg(combinefun)
  
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
    
 # browser()
    # First union or intersection:
    if (combinefun == "unionmodel"){
      mod_union <- unionmodel(mod_prune, matrices = matrices)
    } else if (combinefun == "intersectionmodel"){
      mod_union <- intersectionmodel(mod_prune, matrices = matrices)
    } else {
      
      # for each parameter in the relevant matrices, make equal if both included:
      mod_union <- mod_prune
      
      # But obtain starting values from unionmodel function:
      mod_union@parameters$est <- intersectionmodel(mod_prune, matrices = matrices, runmodel = TRUE) %>%
        groupequal(matrices) %>%
        runmodel %>% 
        '@'('parameters') %>% '$'('est')
      
    }
    
    # # Then set equal:
    # for (m in seq_along(matrices)){
    #   mod_union <- groupequal(mod_union,matrices[m], verbose = FALSE)  
    # }
    # FIXME: Set only the parameters equal that are included in all groups:
    for (m in seq_along(matrices)){
      for (p in which(mod_union@parameters$matrix == matrices[m] & mod_union@parameters$group_id==1)){
        if (!any(mod_union@parameters$fixed[mod_union@parameters$matrix == mod_union@parameters$matrix[p] &
                                            mod_union@parameters$row == mod_union@parameters$row[p] &
                                            mod_union@parameters$col == mod_union@parameters$col[p]])){
          mod_union <- mod_union %>% groupequal(
            matrix = mod_union@parameters$matrix[p],
            row = mod_union@parameters$row[p],
            col = mod_union@parameters$col[p], verbose = FALSE)
        }
        
      }
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
    
    # Final prune step:
    if (final_prune=="partialprune"){
      mod_partialpooled <- curMod %>% prune(alpha=alpha,verbose=FALSE,runmodel=TRUE,matrices=matrices,...)   
    } else {
      
     # loop over all parameters in the matrices:
      all_pars <- which(curMod@parameters$matrix%in% matrices)
      
      # Start loop:
      it <- 1
      while(length(all_pars) >= 1){
        # Parameter:
        p <- all_pars[1]
        
        # If already fixed skio:
        if (curMod@parameters$fixed[p]){
          all_pars <- all_pars[-1]
          next
        }
        
        # If constrained equal, fix all groups to zero if not significant at alpha:
        all_groups_pars <- which(curMod@parameters$matrix ==  curMod@parameters$matrix[p] & 
                                   curMod@parameters$row == curMod@parameters$row[p] & 
                                   curMod@parameters$col == curMod@parameters$col[p])
        
        # Equal constrained?
        eq_constrained <- length(unique(curMod@parameters$par[all_groups_pars])) == 1
        
        # If equal:
        if (eq_constrained){
          # Non-significant?
          sig <- curMod@parameters$p[p] < alpha
          if (!sig){
            curMod <- curMod %>% fixpar(curMod@parameters$matrix[p],curMod@parameters$row[p],curMod@parameters$col[p])
          }
          
          # Remove all parameters from the vector:
          all_pars <- all_pars[!all_pars %in% all_groups_pars]
        } else {
          # If not equal, fix to zero if par was removed in the pruned model:
          if (mod_prune@parameters$fixed[p]){
            curMod <- curMod %>% fixpar(matrix = curMod@parameters$matrix[p],
                                        row = curMod@parameters$row[p],
                                        col = curMod@parameters$col[p],
                                        group = curMod@parameters$group_id[p])
          }
          
          # Remove parameter from the vector:
          all_pars <- all_pars[all_pars != p]
        }
      }
      
      # Run the model again:
      mod_partialpooled <- curMod %>% runmodel(..., verbose = verbose)
    }
    
    
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