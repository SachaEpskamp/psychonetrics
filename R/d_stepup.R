# Stepwise up search:
stepup <- function(
  x, # psychonetrics model
  alpha = 0.01, # Alpha to use for modification indices
  criterion = "bic", # Stop when criterion is no longer improved. Can also be none to ignore
  matrices, # Matrices to search
  mi = c("mi","mi_free","mi_equal"),
  ... # Fit arguments
){
  mi <- match.arg(mi)
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
    } else if (x@model == "lvm"){
      
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
    }  else if (x@model == "var1"){
      matrices <- c("beta")
      if (x@types$zeta == "prec"){
        matrices <- c(matrices,"kappa_zeta")
      } else if (x@types$zeta == "ggm"){
        matrices <- c(matrices,"omega_zeta")
      }
      
    }else if (x@model == "panelvar1"){
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
      
    } else if (x@model == "dlvm1"){
      matrices <- c("beta")
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
      
    }  else stop("No default argument for 'matrices' for current model.")
  }
  
  # Check if MIs are added:
  if (all(is.na(x@parameters[[mi]]))){
    x <- x %>% addMIs(matrices = matrices)
  }
  
  # Start loop:
  repeat{
    oldMod <- x
    # Stepwise up?
    # if (!any(matrices %in% x@equal)){ # FIXME: this will add equality constraints for all matrices...
    if (any(x@parameters[[mi]][x@parameters$matrix %in% matrices & x@parameters$fixed] > qchisq(alpha,1,lower.tail=FALSE))){
      
      # FIXME: Make nice free parameter function
      best <- which(x@parameters[[mi]] == max(x@parameters[[mi]][x@parameters$matrix %in% matrices & x@parameters$fixed]))
      x@parameters$par[best] <- max(x@parameters$par) + 1
      x@parameters$fixed[best] <- FALSE
      
      # Perturb estimate a bit:
      x@parameters$est[best] <- 0.01
      
      # Update the model:
      x@extramatrices$M <- Mmatrix(x@parameters) # FIXME: Make nice function for this
      
      # Run:
      x <- x %>% runmodel(...,log=FALSE)
    } else {
      break
    }
    # } 
    # else {
    #   if (any(x@parameters[[mi]]_equal[x@parameters$matrix %in% matrices & x@parameters$fixed] > qchisq(alpha,1,lower.tail=FALSE))){
    #     # FIXME: Make nice free parameter function
    #     best <- which(x@parameters[[mi]]_equal == max(x@parameters[[mi]]_equal[x@parameters$matrix %in% matrices & x@parameters$fixed]))
    #     x@parameters$par[best] <- max(x@parameters$par) + 1
    #     x@parameters$fixed[best] <- FALSE
    #     
    #     # Perturb estimate a bit:
    #     x@parameters$est[best] <- 0.01
    #     
    #     # Update the model:
    #     x@extramatrices$M <- Mmatrix(x@parameters)
    #     
    #     # Run:
    #     x <- x %>% runmodel(...,log=FALSE)
    #   } else {
    #     break
    #   }
    # }
    # 
    # Check criterion:
    if (criterion != "none"){
      oldCrit <- oldMod@fitmeasures[[criterion]]
      newCrit <- x@fitmeasures[[criterion]]
      if (oldCrit < newCrit){
        x <- oldMod
        break
      }
    }
  }
  
  # Add log:
  x <- addLog(x, "Performed step-up model search")
  
  return(x)
  
}