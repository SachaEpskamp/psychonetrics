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
  useMIs = c("joint","simple"),
  release_model = c("pruned","saturated"),
  identify = TRUE,
  ...){
  useMIs <- match.arg(useMIs)
  release_model <- match.arg(release_model)
  # Verbose:
  if (missing(verbose)){
    verbose <- x@verbose
  }
  fixed <- NULL
  mi_free <- NULL
  equal <- NULL
  
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
      
    } else if (x@model == "panelvar"){
      matrices <- c("beta")
      if (x@types$within_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_within")
      } else if (x@types$within_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_within")
      }

      if (x@types$between_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_between")
      } else if (x@types$between_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_between")
      }

    } else if (x@model == "ml_varcov"){
      matrices <- character(0)
      if (x@types$within == "prec"){
        matrices <- c(matrices,"kappa_within")
      } else if (x@types$within == "ggm"){
        matrices <- c(matrices,"omega_within")
      } else if (x@types$within == "cor"){
        matrices <- c(matrices,"rho_within")
      } else if (x@types$within == "chol"){
        matrices <- c(matrices,"lowertri_within")
      } else {
        matrices <- c(matrices,"sigma_within")
      }
      if (x@types$between == "prec"){
        matrices <- c(matrices,"kappa_between")
      } else if (x@types$between == "ggm"){
        matrices <- c(matrices,"omega_between")
      } else if (x@types$between == "cor"){
        matrices <- c(matrices,"rho_between")
      } else if (x@types$between == "chol"){
        matrices <- c(matrices,"lowertri_between")
      } else {
        matrices <- c(matrices,"sigma_between")
      }
    } else if (x@model == "ml_var1"){
      matrices <- c("beta")
      if (x@types$within_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_within")
      } else if (x@types$within_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_within")
      }

      if (x@types$between_latent == "prec"){
        matrices <- c(matrices,"kappa_zeta_between")
      } else if (x@types$between_latent == "ggm"){
        matrices <- c(matrices,"omega_zeta_between")
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
      
    }  else if (x@model %in% c("Ising", "BlumeCapel")){
      matrices <- c("omega")
      
    }  else stop("No default argument for 'matrices' for current model.")

    # PDC temporal parameterization: the temporal matrix is "PDC", not "beta":
    if ((!is.null(x@types$temporal) && x@types$temporal == "PDC") ||
        (!is.null(x@types$temporal_latent) && x@types$temporal_latent == "PDC")){
      matrices[matrices == "beta"] <- "PDC"
    }
  }
  
  # Prune first:
  # NOTE: the trailing runmodel looks redundant (prune(runmodel = TRUE) already
  # runs the pruned model), but it re-polishes the optimum and removing it
  # changes estimates at ~1e-8, so it is kept for exact backward compatibility:
  if (verbose) message("Pruning model...")
  mod_prune <- prune(x,alpha=alpha,verbose=FALSE,runmodel=TRUE,matrices=matrices,identify=identify,...) %>% runmodel
  
  # if not empty, look for equality:
  if (!all(mod_prune@parameters$est[mod_prune@parameters$matrix%in%matrices&!mod_prune@parameters$fixed]==0)){
    
    
    # Combine models:
    if (verbose) message("Combining models...")
    
 # browser()
    # First union or intersection:
    if (combinefun == "unionmodel"){
      mod_union <- unionmodel(mod_prune, matrices = matrices,identify =identify)
    } else if (combinefun == "intersectionmodel"){
      mod_union <- intersectionmodel(mod_prune, matrices = matrices,identify =identify)
    } else {
      
      # for each parameter in the relevant matrices, make equal if both included:
      mod_union <- mod_prune
      
      # But obtain starting values from unionmodel function:
      mod_union@parameters$est <- intersectionmodel(mod_prune, matrices = matrices, runmodel = TRUE,identify =identify ) %>%
        groupequal(matrices, identify = identify) %>%
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
            col = mod_union@parameters$col[p], verbose = FALSE,identify =identify )
        }
        
      }
    }
  
    # Then run:
    mod_union <- runmodel(mod_union, verbose = FALSE)
    
    
    
    if (verbose) message("Partial pruning...")
    curMod <- mod_union
    
    
    repeat{
      pars <- curMod@parameters
      pars$equal <- pars$par != 0 & (duplicated(pars$par) | rev(duplicated(rev(pars$par))))

      if (useMIs == "joint"){
        # Use the joint score-test (Lagrange multiplier) statistic, one value per
        # equality-constrained (matrix,row,col). Stored on every group's row; take the
        # first to avoid double counting.
        miDF <- pars %>% filter(drop(!fixed), drop(equal)) %>%
          filter(drop(matrix %in% matrices)) %>%
          group_by(.data[["row"]],.data[["col"]],.data[["matrix"]]) %>%
          summarize(mi_free = dplyr::first(.data[["mi_free_joint"]])) %>%
          filter(drop(!is.na(.data[["mi_free"]]))) %>%
          arrange(-mi_free)
      } else {
        # Legacy (<= 0.15.3): sum of per-group univariate score statistics.
        miDF <- pars %>% filter(drop(!fixed), drop(equal)) %>% group_by(.data[["row"]],.data[["col"]],.data[["matrix"]]) %>%
          filter(drop(matrix %in% matrices)) %>%
          summarize(mi_free = sum(mi_free)) %>%
          arrange(-mi_free)
      }
      
      # break if empty:
      if (nrow(miDF)==0){
        break
      }


            
      # Free the best parameter:
      tryres <- try({
      mat_i <- miDF$matrix[1]; row_i <- miDF$row[1]; col_i <- miDF$col[1]
      propMod <- curMod %>%
        groupfree(mat_i, row_i, col_i,identify =identify )
      # release_model = "pruned": before refitting and comparing BIC, fix the
      # just-released parameter to zero in any group where Step 1's mod_prune
      # had it removed. This tests the equality constraint against the
      # asymmetric structure suggested by the per-group prune (rather than
      # against the fully saturated per-group alternative), which improves
      # specificity at low sample sizes.
      if (release_model == "pruned"){
        prune_rows <- which(mod_prune@parameters$matrix == mat_i &
                              mod_prune@parameters$row == row_i &
                              mod_prune@parameters$col == col_i)
        for (pr in prune_rows){
          if (isTRUE(mod_prune@parameters$fixed[pr])){
            propMod <- propMod %>% fixpar(
              matrix = mat_i, row = row_i, col = col_i,
              group = mod_prune@parameters$group_id[pr],
              verbose = FALSE)
          }
        }
      }
      # Inner-loop refit: only need BIC (from addfit) and the joint score-test
      # statistic mi_free_joint to pick the next release. Skip addMIs (which
      # calls psychonetrics_FisherInformation_cpp three times for type =
      # normal/free/equal). KEEP addSEs and addInformation enabled because
      # the post-loop prune() reads p-values (from addSEs) and @information
      # (from addInformation) on the final curMod. Disabling either breaks
      # the post-loop prune step on cases where the inner loop iterates
      # more than once.
      # For multi-group Ising at high node counts the addMIs Fisher builds
      # dominate the per-iteration cost because each Fisher build enumerates
      # 2^N states.
      propMod <- propMod %>% runmodel(addMIs = FALSE)
      # Re-populate mi_free_joint via the joint score test only (the inner
      # loop reads this column to choose which equality constraint to release
      # next). Mirrors the joint-test write block in addMIs_inner_full(type =
      # "free"). For useMIs == "simple" we fall back to a full addMIs call
      # to populate the legacy per-parameter mi_free column the simple
      # branch reads.
      if (useMIs == "joint" && nrow(propMod@sample@groups) >= 2){
        if (is.null(propMod@parameters$mi_free_joint))
          propMod@parameters$mi_free_joint <- NA_real_
        if (is.null(propMod@parameters$pmi_free_joint))
          propMod@parameters$pmi_free_joint <- NA_real_
        if (is.null(propMod@parameters$df_free_joint))
          propMod@parameters$df_free_joint <- NA_real_
        propMod@parameters$mi_free_joint[]  <- NA_real_
        propMod@parameters$pmi_free_joint[] <- NA_real_
        propMod@parameters$df_free_joint[]  <- NA_real_
        est_res <- tryCatch(
          .equalityScoreTestInner(propMod, method = "jacobian",
                                  joint_only = TRUE),
          error = function(e) NULL
        )
        if (!is.null(est_res) && !is.null(est_res$total) && nrow(est_res$total) > 0){
          tot <- est_res$total
          for (ii in seq_len(nrow(tot))){
            wr <- which(propMod@parameters$matrix == tot$matrix[ii] &
                          propMod@parameters$row == tot$row[ii] &
                          propMod@parameters$col == tot$col[ii])
            propMod@parameters$mi_free_joint[wr]  <- round(tot$X2[ii], 10)
            propMod@parameters$pmi_free_joint[wr] <- round(tot$p.value[ii], 10)
            propMod@parameters$df_free_joint[wr]  <- tot$df[ii]
          }
        }
      } else if (useMIs == "simple"){
        propMod <- addMIs(propMod, verbose = FALSE)
      }
      })
      
      # if (inherits(tryres, "try-error")){
      #   browser()
      # }
      
      # Test BIC:
      if (propMod@fitmeasures$objective > 0 && propMod@fitmeasures$bic < curMod@fitmeasures$bic){
        curMod <- propMod
      } else {
        break
      }
    }
    
    # Final prune step:
    if (final_prune=="partialprune"){
      mod_partialpooled <- curMod %>% prune(alpha=alpha,verbose=FALSE,runmodel=TRUE,matrices=matrices,identify=identify,...)
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
          if (!is.na(sig) && !sig){
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