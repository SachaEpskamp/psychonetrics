# Function to prune all non-significant results:
unionmodel <- function(
  x, # Model
  runmodel = FALSE,
  verbose,
  log = TRUE,
  identify = TRUE,
  matrices, # Automatically chosen
  ...){
  
  stopifnot(is(x,"psychonetrics"))
  
  # Choose matrices:
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
  
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  # If only one group, nothing to do!
  if (nrow(x@sample@groups) == 1){
    if (verbose) message("Only one group in model, unionmodel does nothing...")
    return(x)
  }
  
  # Copy parameter table:
  pars <- x@parameters
  
  # Add dummy id:
  pars$id <- seq_len(nrow(pars))
  
  # Find set of parameters at least in one group:
  pars <- pars %>%  left_join(pars %>% dplyr::group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% dplyr::summarise(anyFree = any(!.data[['fixed']])),
                    by = c("matrix","row","col"))
  
  # So which to free?
  whichFree <- which(pars$fixed & pars$anyFree & pars$matrix %in% matrices)
  
  # If nothing to do, return
  if (length(whichFree)==0){
    return(x)
  }
  

  # Set computed:
  x@computed <- FALSE
  
  # Free the parameters:
  x@parameters$fixed[whichFree] <- FALSE
  x@parameters$par[whichFree] <- max(x@parameters$par) + seq_len(length(whichFree))
  # x@parameters$mi[whichFree] <- NA
  # x@parameters$pmi[whichFree] <- NA
  # x@parameters$mi_equal[whichFree] <- NA
  # x@parameters$pmi_equal[whichFree] <- NA
  
  x@parameters <- clearpars(x@parameters,whichFree)

  # Just in case... # FIXME: Not needed right?
  x@parameters   <- parRelabel(x@parameters)
  
  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  if (verbose){
    message(paste0("Freeing ",length(whichFree)," parameters!"))
  }
  
  
  if (log){
    # Add log:
    x <- addLog(x, paste("Unified models across groups (union model). Freed in total ",length(whichFree)," parameters")) 
  }

  # Rerun if needed:
  if (runmodel){
    x <- x %>% runmodel(verbose=verbose,...)
  }

  x
}