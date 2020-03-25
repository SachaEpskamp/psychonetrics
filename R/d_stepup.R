# Stepwise up search:
stepup <- function(
  x, # psychonetrics model
  alpha = 0.01, # Alpha to use for modification indices
  criterion = "bic", # Stop when criterion is no longer improved. Can also be none to ignore
  matrices, # Matrices to search
  mi = c("mi","mi_free","mi_equal"),
  greedyadjust = c("bonferroni", "none", "holm", "hochberg", "hommel", "fdr", "BH", "BY"),
  stopif,
  greedy = FALSE, # If TRUE, will start by adding all significant effects followed by pruning
  verbose,
  checkinformation = TRUE,
  singularinformation = c("tryfix","skip","continue","stop"), # tryfix = try to fix by adjusting starting values (once), skip = go to next parameter, continue = continue search, stop = stop and return current and previous model
  startEPC = TRUE,
    # maxtry = 0,
  ... # Fit arguments
){
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  maxtry <- 1
  singularinformation <- match.arg(singularinformation)
  triedfixing <- FALSE
  
  # If not run, run model:
  if (!x@computed){
    x <- x %>% runmodel(..., verbose = verbose)
  }
  
  greedyadjust <- match.arg(greedyadjust)
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
  
  # Check if MIs are added:
  if (all(is.na(x@parameters[[mi]]))){
    x <- x %>% addMIs(matrices = matrices)
  }
  
  # Start loop:
  repeat{
    oldMod <- x
    # Stepwise up?
    # if (!any(matrices %in% x@equal)){ # FIXME: this will add equality constraints for all matrices...
    if (any(x@parameters[[mi]][x@parameters$matrix %in% matrices & x@parameters$fixed & !x@parameters$identified & !is.na(x@parameters[[mi]])] > qchisq(alpha,1,lower.tail=FALSE))){
      
      # FIXME: Make nice free parameter function
      if (!greedy){
        x@parameters[[mi]] <- ifelse(is.na(x@parameters[[mi]]),0,x@parameters[[mi]])
        best <- which(x@parameters$matrix %in% matrices & x@parameters[[mi]] == max(x@parameters[[mi]][x@parameters$matrix %in% matrices & x@parameters$fixed & !x@parameters$identified]))[1]
        
        # Check if equal constrained:
        if (x@parameters$matrix[best] %in% x@equal){
          x <- freepar(x, matrix = x@parameters$matrix[best],row = x@parameters$row[best],
                       col = x@parameters$col[best], 
                       verbose = FALSE, log = FALSE, startEPC=startEPC)
          x <- groupequal(x, matrix = x@parameters$matrix[best],row = x@parameters$row[best],
                          col = x@parameters$col[best],verbose = FALSE, log = FALSE)
          
          if (verbose){
            message(paste0("Adding parameter ",x@parameters$matrix[best],"[",x@parameters$var1[best],", ",x@parameters$var2[best],"]"))              
          }
          
          
        } else {
          x <- freepar(x, matrix = x@parameters$matrix[best],row = x@parameters$row[best],
                       col = x@parameters$col[best], group = x@parameters$group_id[best],
                       verbose = FALSE, log = FALSE, startEPC=startEPC)
          
          if (verbose){
            if (nrow(x@sample@groups) == 1){
              message(paste0("Adding parameter ",x@parameters$matrix[best],"[",x@parameters$var1[best],", ",x@parameters$var2[best],"]"))            
            } else {
              message(paste0("Adding parameter ",x@parameters$matrix[best],"[",x@parameters$var1[best],", ",x@parameters$var2[best],", ",x@parameters$group[best],"]"))
            }
            
          }
        }
        
        # x@parameters$par[best] <- max(x@parameters$par) + 1
        # x@parameters$fixed[best] <- FALSE
        # 
        # # # Perturb estimate a bit:
        # # x@parameters$est[best] <- 0.01
        # # Set estimate to EPC:
        # x@parameters$est[best] <- x@parameters$est[best] + sign(x@parameters$epc[best]) * 0.01
        # 
        # # Update the model:
        # x@extramatrices$M <- Mmatrix(x@parameters) # FIXME: Make nice function for this
        # 
        
        
        # Run:
        curtry <- 0
        repeat{
          suppressWarnings(newx <- x %>% runmodel(...,log=FALSE))
          
          # Check information:
          if (checkinformation){
            # if (any(eigen(newx@information)$values < -sqrt(.Machine$double.eps))){
            if (!sympd_cpp(newx@information)){
              # if (curtry < maxtry){
              #   if (verbose){
              #     message(paste("Model may not be identified, adjusting start values and trying again."))
              #   }
              #   x@parameters$est[best] <- 0.1 * x@parameters$est[best]
              #   curtry <- curtry + 1
              # } else {
              # if (verbose){
              #   message(paste("Model may not be identified, continuing with previous model."))
              # }
              # x <- oldMod
              # x@parameters$identified[best] <- TRUE
              # x@parameters[[mi]][best] <- 0
              # break
              # }
              
              if (singularinformation == "tryfix"){
                if (triedfixing){
                  message("Could not repair identification issue. Aborting search and returning previous model.")
                  return(oldMod)
                }
                
                triedfixing <- TRUE
                
                if (verbose){
                  message("Model may not be identified. Adjusting starting values and trying again")
                }
                x <- runmodel(emergencystart(x))
                
              } else if (singularinformation == "skip"){
                if (verbose){
                  message(paste("Model may not be identified, continuing with previous model."))
                }
                x <- oldMod
                x@parameters$identified[best] <- TRUE
                x@parameters[[mi]][best] <- 0
                break
              } else if (singularinformation == "stop"){
                if (verbose){
                  message(paste("Model may not be identified, returning list with last two models."))
                }
                return(list(
                  final = newx,
                  previous = x
                ))
              }
              
              
              
            } else {
              triedfixing <- FALSE
              x <- newx
              break
            }
          } else {
            x <- newx
            break
          }
        }
        
        # Check if fit actually improved:
        compare <- compare(old = oldMod, new = x)
        if (!compare$p_value[2] < alpha){
          if (verbose){
            message(paste("Model did not improve at given alpha, returning previous model."))
          }
          x <- oldMod
          break
        }
        
        # Check criterion:
        if (criterion != "none"){
          if (!criterion %in% names(oldMod@fitmeasures)){
            stop(paste0("Criterion '",criterion,"' is not supported."))
          }
          oldCrit <- oldMod@fitmeasures[[criterion]]
          newCrit <- x@fitmeasures[[criterion]]
          
          if (oldCrit < newCrit){
            if (verbose){
              message(paste("Model did not improve criterion, returning previous model."))
            }
            x <- oldMod
            break
          }
        }
        
      } else {
        # Add al significant effects:
        # nTest <- sum(x@parameters$matrix %in% matrices & x@parameters$fixed)
        # best <- which(x@parameters[[mi]] > qchisq(alpha,1,lower.tail=FALSE) & x@parameters$matrix %in% matrices & x@parameters$fixed)
        
        parsToTest <- which(x@parameters$matrix %in% matrices & x@parameters$fixed & !x@parameters$identified)
        x@parameters[[mi]] <- ifelse(is.na(x@parameters[[mi]]),0,x@parameters[[mi]])
        best <- parsToTest[p.adjust(pchisq(x@parameters[[mi]][parsToTest],1,lower.tail=FALSE), method = greedyadjust) < alpha]
        
        x@parameters$par[best] <- max(x@parameters$par) + seq_along(best)
        x@parameters$fixed[best] <- FALSE
        
        # Perturb estimate a bit:
        # x@parameters$est[best] <- 0.001
        # Set estimate to EPC:
        x@parameters$est[best] <- x@parameters$est[best] + sign(x@parameters$epc[best]) * 0.01
        # Use the emergency start function:
        x <- emergencystart(x)
        
        # Update the model:
        # x@extramatrices$M <- Mmatrix(x@parameters) # FIXME: Make nice function for this
        if (x@cpp){
          x@extramatrices$M  <- Mmatrix_cpp(x@parameters)
        } else {
          x@extramatrices$M  <- Mmatrix(x@parameters)  
        }
        
        if (verbose){
          message(paste("Adding",length(best),"parameters in greedy search start."))
        }
        # Run:
        greedy <- FALSE
        curtry <- 0
        repeat{
          newx <- x %>% runmodel(...,log=FALSE) # %>% prune(alpha = alpha, adjust = greedyadjust)
          
          if (checkinformation){
            
            # if (any(eigen(x@information)$values < -sqrt(.Machine$double.eps))){
            if (!sympd_cpp(x@information)){
              if (curtry < maxtry){
                if (verbose){
                  message(paste("Model may not be identified, adjusting start values and trying again."))
                }
                x <- emergencystart(x)
                curtry <- curtry + 1
              } else {
                if (verbose){
                  message(paste("Model may not be identified, continuing with previous model."))
                }
                x <- oldMod  
                break
              }
            } else {
              x <- newx #%>% prune(alpha = alpha, adjust = greedyadjust)
              break
            }
          } else {
            x <- newx #%>% prune(alpha = alpha, adjust = greedyadjust)
            break
          }
          
        }
        
        # Check if fit actually improved:
        compare <- compare(old = oldMod, new = x)
        if (!compare$p_value[2] < alpha){
          if (verbose){
            message(paste("Model did not improve at given alpha, returning previous model."))
          }
          x <- oldMod
          break
        }
      
      
        # Check criterion:
        if (criterion != "none"){
          if (!criterion %in% names(oldMod@fitmeasures)){
            stop(paste0("Criterion '",criterion,"' is not supported."))
          }
          oldCrit <- oldMod@fitmeasures[[criterion]]
          newCrit <- x@fitmeasures[[criterion]]
          
          if (oldCrit < newCrit){
            if (verbose){
              message(paste("Model did not improve criterion, continuing with previous model."))
            }
            x <- oldMod
          }
        }
        
      }
      
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
    
    # Check if we should stop
    if (!missing(stopif)){
      condition <- eval(substitute(stopif), envir = x@fitmeasures)
      if (condition){
        if (verbose){
          message(paste("Model in line with stopping criterion: stopping search."))
        }
        break
      }
    }
    
  }
  
  # Add log:
  x <- addLog(x, "Performed step-up model search")
  
  return(x)
  
}