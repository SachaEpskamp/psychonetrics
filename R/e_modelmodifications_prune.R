# Function to prune all non-significant results:
prune <- function(
  x, # Model
  alpha = 0.01, # Significance
  # bonferroni = FALSE,
  adjust = c( "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
  matrices, # Automatically chosen
  runmodel = TRUE,
  recursive = FALSE,
  verbose,
  log = TRUE,
  identify = TRUE,
  # nCores = 1,
  # reps = 1000,
  startreduce = 1,
  limit = Inf,
  ...){
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  reps <- 1000
  nCores <- 1
  bootstrap <- FALSE
  adjust <- match.arg(adjust)
  # If not run, run model:
  if (!x@computed){
    x <- x %>% runmodel(..., verbose = verbose)
  }
  
  if (bootstrap && all(is.na(x@parameters$boot_p))){
    if (verbose) message("Bootstrapping SEs...")
    x <- x %>% bootstrap_SEs(nCores = nCores, reps = reps,verbose = verbose)
  }
  
  # Whcih cols?
  if (bootstrap){
    secol <- "se_boot"
    pcol <- "p_boot" 
  } else {
    secol <- "se"
    pcol <- "p"
  }
  
  
  if (!runmodel & recursive){
    stop("recursive = TRUE requires runmodel = TRUE")
  }
  # Select the matrices:
  # if (missing(matrices)){
  #   if (x@model == "ggm"){
  #     matrices <- "omega"
  #   } else if (x@model == "precision"){
  #     matrices <- "kappa"
  #   } else if (x@model == "lnm"){
  #     matrices <- "omega_eta"
  #   } else if (x@model == "rnm"){
  #     matrices <- "omega_epsilon"
  #   } else if (x@model == "gvar"){
  #     matrices <- c("beta","omega_zeta")
  #   } else if (x@model == "varcov"){
  #     matrices <- "sigma"
  #   } else if (x@model == "cholesky"){
  #     matrices <- "lowertri"
  #   } else {
  #     stop("'no default 'matrix' argument implemented yet..")
  #   }
  # }
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
  
  # Which parameters to test:
  whichTest <- which(x@parameters$matrix %in% matrices & !x@parameters$fixed & (x@parameters$var1_id!=x@parameters$var2_id | x@parameters$matrix == "beta"))
  
  # Number of tests:
  nTest <- length(unique(x@parameters$par[whichTest]))
  
  # If no tests, break:
  if (nTest == 0){
    return(x)
  }
  
  # Adjust alpha:
  # if (bonferroni){
  #   alpha_adjust <- alpha / nTest    
  # } else {
  #   alpha_adjust <- alpha
  # }
  # Test for significance:
  # nonsig <- x@parameters$p > alpha_adjust & (seq_len(nrow(x@parameters)) %in% whichTest)
  pValues <- p.adjust(x@parameters[[pcol]],method = adjust) 
  nonsig <- p.adjust(x@parameters[[pcol]],method = adjust) > alpha & (seq_len(nrow(x@parameters)) %in% whichTest)
  
  
  # If any non-sig, adjust:
  if (all(is.na(nonsig)) || !any(nonsig[!is.na(nonsig)])){
    if (log){
      # Add log:
      x <- addLog(x, paste("Pruned all parameters in matrices ",paste0(matrices,collapse=", ")," at alpha = ",alpha," (none were removed)")) 
    }    
    return(x)
  } 
  
  # Limit:
  nonsig[nonsig & pValues < min(head(sort(pValues[nonsig],decreasing = TRUE), limit))] <- FALSE
  
  
  curPars <- max(x@parameters$par)
  
  # Set computed:
  x@computed <- FALSE
  
  # Clear the parameters:
  x@parameters$est[nonsig] <- 0
  # x@parameters$std[nonsig] <- 0
  # x@parameters$se[nonsig] <- 0
  x@parameters$fixed[nonsig] <- TRUE
  x@parameters$par[nonsig] <- 0
  # x@parameters$mi[nonsig] <- NA
  # x@parameters$pmi[nonsig] <- NA
  # x@parameters$mi_equal[nonsig] <- NA
  # x@parameters$pmi_equal[nonsig] <- NA
  
  # FIXME: Reduce parameter estimates from remainder of matrix a bit to avoid problems:
  x@parameters$est[x@parameters$matrix %in% matrices & !x@parameters$fixed & !x@parameters$identified & x@parameters$est != 0] <- 
    startreduce * x@parameters$est[x@parameters$matrix %in% matrices & !x@parameters$fixed & !x@parameters$identified & x@parameters$est != 0] 
  
  x@parameters <- clearpars(x@parameters,nonsig)
  
  x@parameters   <- parRelabel(x@parameters)
  
  # Identify:
  if (identify){
    x <- identify(x)
  }
  
  if (verbose){
    message(paste0("Clearing ",curPars - max(x@parameters$par)," parameters!"))
  }
  
  
  if (log){
    # Add log:
    # x <- addLog(x, paste("Pruned all parameters in matrices ",paste0(matrices,collapse=", ")," at alpha = ",alpha_adjust))
    x <- addLog(x, paste("Pruned all parameters in matrices ",paste0(matrices,collapse=", ")," at alpha = ",alpha,ifelse(adjust=="none","",paste0(" (adjusted: ",adjust,")"))) )
  }
  
  # Rerun if needed:
  if (runmodel){
    xOld <- x
    suppressWarnings(x <- x %>% runmodel(verbose=verbose,...))
    
    # If not identified, try with emergency start:
    # if (any(eigen(x@information)$values < -sqrt(.Machine$double.eps))){
    if (!sympd_cpp(x@information)){
      # cat("EMERGENCYSTART")
      x <- emergencystart(xOld) %>% runmodel(verbose=verbose,...)
    }
  }
  
  # Recurse if needed:
  if (recursive){
    x <- x %>% prune(
      alpha = 0.01, # Significance
      adjust = adjust,
      matrices = matrices, # Automatically chosen
      runmodel = TRUE,
      recursive = TRUE,
      verbose = verbose,
      log = log,
      identify = identify,
      # bootstrap = bootstrap,
      ...
    )
  }
  
  
  x
}