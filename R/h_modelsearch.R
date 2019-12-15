# ggmModSelect inspired algorithm:
modelsearch <- function(x,
                        criterion = "bic", # Stop when criterion is no longer improved. Can also be none to ignore
                        matrices, # Matrices to search
                        prunealpha = 0.001, # Minimum p-value for edges tested to be removed
                        addalpha = 0.05, # Maximum p-value for edges tested to be added
                        verbose = TRUE,
                        ...
){
  mi <- "mi"
  pmi <- "pmi"
  
  # FIXME: If number of groups > 1, stop:
  if (nrow(x@sample@groups) > 1) stop("'modelsearch' is only implemented for single group models at the moment.")
  
  # If not run, run model:
  if (!x@computed){
    x <- x %>% runmodel(..., verbose = verbose)
  }
  
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
    } else if (x@model == "meta_varcov"){
      
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
      
    }  else stop("No default argument for 'matrices' for current model.")
  }
  
  # Check if MIs are added:
  if (all(is.na(x@parameters[[mi]]))){
    x <- x %>% addMIs(matrices = matrices)
  }
  
  # Start outer loop for evaluating ALL competing models!
  repeat{
    
    # List all parameters in matrices of interest:
    ind <- which(x@parameters$matrix %in% matrices & !x@parameters$identified)
    
    allParsToConsider <- ind[
      # Non significant edges to remove?
      (!x@parameters$fixed[ind]  & x@parameters$p[ind] > prunealpha) |
        
        # Significant edges to add?
        (x@parameters$fixed[ind]  & x@parameters[[pmi]][ind] < prunealpha)
      ]
    
 
    if (length(allParsToConsider) == 0){
      if (verbose){
        message("No more models to search at given alpha levels.")
        return(x)
      }
    }
    
    # Set parameters to consider:
    parsToConsider <- allParsToConsider
    
    # Start the inner loop for subset of models:
    repeat{
      
      # Set progress bar:
      if (verbose){
        if (length(allParsToConsider) == length(parsToConsider)){
          message("Evaluating all models at given alpha levels...")
        } else {
          message("Evaluating subset of models...")
        }
        
        pb <- txtProgressBar(min = 0, max = length(parsToConsider), initial = 0, style = 3)
      }
      
      # Empty list of models:
      propMods <- vector("list", length(parsToConsider))
      
      # Loop over these parameters:
      for (i in seq_along(parsToConsider)){
        curpar <- parsToConsider[i]
        
        # Existing parameter?
        fixed <- x@parameters$fixed[curpar]
        
        if (fixed){
          # Try to add (stepup)
          propMods[[i]] <- x %>% freepar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                         group = x@parameters$group_id[curpar],
                                         start = x@parameters$epc[curpar], verbose = FALSE)
          
          # Run the model:
          propMods[[i]] <- propMods[[i]] %>% runmodel(verbose = FALSE, addMIs = FALSE, ...)
          
          # Fisher information ok?
          ev <- eigen(propMods[[i]]@information)$values
          
          # If not, try again with different starts:
          if (any(ev < -sqrt(.Machine$double.eps))){
            try <- 1
            repeat{
              if (try == 1){
                # First try, 0.5 * EPC:
                propMods[[i]] <- x %>% freepar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                               group = x@parameters$group_id[curpar],
                                               start = 0.5 * x@parameters$epc[curpar], verbose = FALSE)
                
              } else if (try == 2){
                # Second try, 0.0001 * sign(EPC):
                propMods[[i]] <- x %>% freepar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                               group = x@parameters$group_id[curpar],
                                               start = 0.0001 * sign(x@parameters$epc[curpar]), verbose = FALSE)
                
              } else if (try == 3){
                # Final try, emergency start:
                propMods[[i]] <- x %>% freepar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                               group = x@parameters$group_id[curpar],
                                               start = x@parameters$epc[curpar], verbose = FALSE) %>% emergencystart
                
              } else if (try > 3){
                
                # Give up:
                propMods[[i]] <- x
                propMods[[i]]@fitmeasures[[criterion]] <- Inf
                break
                
              }
              
              # Run the model:
              propMods[[i]] <- propMods[[i]] %>% runmodel(verbose = FALSE, addMIs = FALSE, ...)
              
              # Fisher information ok?
              ev <- eigen(propMods[[i]]@information)$values
              
              if (all(ev > -sqrt(.Machine$double.eps))){
                break
              } else {
                try <- try + 1
              }
              
              
            }
            
          }
          
          
          
        } else {
          
          # Try to remove (prune):
          propMods[[i]] <- x %>% fixpar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                        group = x@parameters$group_id[curpar], verbose = FALSE)
          
          # Run the model:
          propMods[[i]] <- propMods[[i]] %>% runmodel(verbose = FALSE, addMIs = FALSE, ...)
          
          # Fisher information ok?
          ev <- eigen(propMods[[i]]@information)$values
          
          # If not, try again with different starts:
          if (any(ev < -sqrt(.Machine$double.eps))){
            try <- 1
            repeat{
              if (try == 1){
                # First try, reduce all parameters in matrix a bit:
                propMods[[i]] <- x %>% fixpar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                              group = x@parameters$group_id[curpar], verbose = FALSE)
                
                ind <- propMods[[i]]@parameters$matrix %in% propMods[[i]]@parameters$matrices & !propMods[[i]]@parameters$fixed & !propMods[[i]]@parameters$identified &
                  propMods[[i]]@parameters$row != propMods[[i]]@parameters$col
                
                propMods[[i]]@parameters$est[ind] <- 0.5 * propMods[[i]]@parameters$est[ind]
                
                
              } else if (try == 2){
                # Second try, emergency start:
                propMods[[i]] <- x %>% fixpar(matrix = x@parameters$matrix[curpar], row = x@parameters$row[curpar], col = x@parameters$col[curpar], 
                                              group = x@parameters$group_id[curpar], verbose = FALSE) %>% emergencystart
                
              } else if (try > 2){
                
                # Give up:
                propMods[[i]] <- x
                propMods[[i]]@fitmeasures[[criterion]] <- Inf
                break
                
              }
              
              
              # Run the model:
              propMods[[i]] <- propMods[[i]] %>% runmodel(verbose = FALSE, addMIs = FALSE, ...)
              
              # Fisher information ok?
              ev <- eigen(propMods[[i]]@information)$values
              
              if (all(ev > -sqrt(.Machine$double.eps))){
                break
              } else {
                try <- try + 1
              }
              
              
            }
            
          }
          
          
          
        }
        
        # Update progress:
        if (verbose){
          setTxtProgressBar(pb, i)
        }
        
      }
      
      # Update progress:
      if (verbose){
        close(pb)
      }
      
      # Current criterion:
      curCrit <- x@fitmeasures[[criterion]]
      
      # Criterions of proposal models:
      propCrits <- sapply(propMods, function(xx) xx@fitmeasures[[criterion]] )
      
      # which improve?
      whichImprove <- which(propCrits < curCrit)
      
      # if none, break:
      if (length(whichImprove) == 0){
        break # FIXME
      }
      
      # Which is the best?
      whichBest <- which(propCrits == min(propCrits))
      if (length(whichBest) > 1){
        warning("Multiple equivalent models found, selecting one at random.")
        whichBest <- sample(whichBest,1)
      }
      
      # Say something:
      if (verbose){
        if (propMods[[whichBest]]@fitmeasures$df < x@fitmeasures$df){
          message("Adding one parameter...")
        } else {
          message("Fixing one parameter...")
        }        
      }

      
      # Update the model:
      x <- propMods[[whichBest]] %>% addMIs(verbose = FALSE)
      
      
      # Update list of candidates:
      parsToConsider <- parsToConsider[whichImprove[whichImprove != whichBest]]
      
      # if none, break:
      if (length(parsToConsider) == 0){
        break # FIXME
      }
      
    }
    
    # Did we test all models in the last run? If so, we are done (jaj)!
    if (length(allParsToConsider) == length(parsToConsider)){
      if (verbose) message("No more model found to improve fit. Returning optimal model.")
      break
    }
    
  }
  
  # Return the final model:
  return(x)
}