partialprune <- function(
  x, # Model
  alpha = 0.01, # Significance
  matrices, # Automatically chosen
  verbose,
  combinefun = unionmodel,
  ...){
  # Verbose:
  if (missing(verbose)){
    verbose <- x@verbose
  }
  fixed <- NULL
  mi_free <- NULL
  
  
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
  mod_prune <- prune(x,alpha=alpha,verbose=FALSE,runmodel=TRUE,...) 
  
  
  # Combine models:
  if (verbose) message("Combining models...")
  mod_union <- runmodel(groupequal(combinefun(mod_prune),matrices, verbose = FALSE), verbose = FALSE)
  
  if (verbose) message("Partial pruning...")
  curMod <- mod_union
  
  repeat{
    pars <- curMod@parameters

    # Make a data frame in which the equality free parameters are summed:
    # miDF <- pars %>% filter(!fixed) %>% group_by_("group","row","col","matrix") %>%
    #   summarize(mi_free = sum(mi_free)) %>% 
    #   arrange(-mi_free)
    miDF <- pars %>% filter(!fixed) %>% group_by_("row","col","matrix") %>%
      summarize(mi_free = sum(mi_free)) %>% 
      arrange(-mi_free)
    
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
  mod_partialpooled <- curMod %>% prune(alpha=alpha,verbose=FALSE,runmodel=TRUE,...) 
  
  # Select best model:
  mods <- list(x, mod_prune, mod_union, mod_partialpooled)
  
  best <- which.min(c(
    x@fitmeasures$bic,
    mod_prune@fitmeasures$bic,
    mod_union@fitmeasures$bic,
    mod_partialpooled@fitmeasures$bic
  ))
  
  bestmod <- mods[[best]]
  
  return(bestmod)
}