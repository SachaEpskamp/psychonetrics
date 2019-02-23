# this function will automatically identify models:
identify <- function(x){
  if (x@model == "ggm" | x@model == "precision"){
    # Nothing to do..
    return(x)
  }
  
  if (x@model == "lnm"){
    x <- identify_lnm(x)
    return(x)
  }
  
  stop("Model not supported...")
}

# LNM identifier:
identify_lnm <- function(x){
  # Idenitfy type:
  type <- x@identification
  # Single group is easy:
  if (nrow(x@sample@groups) == 1){
    
    # Set all means to zero:
    means <- which(x@parameters$matrix %in% c("mu_eta"))
    x@parameters$est[means] <- 0
    # x@parameters$std[means] <- NA
    x@parameters$par[means] <- 0
    # x@parameters$se[means] <- NA
    # x@parameters$p[means] <- NA
    # x@parameters$mi[means] <- NA
    # x@parameters$pmi[means] <- NA
    # x@parameters$mi_equal[means] <- NA
    # x@parameters$pmi_equal[means] <- NA
    x@parameters$fixed[means] <- TRUE
    x@parameters$identified[means] <- TRUE
    
    x@parameters <- clearpars(x@parameters, means)
    
    
    # Scaling ifentification:
    if (type == "scaling"){
      # Set al scaling to 1:
      scaling <- which(x@parameters$matrix %in% c("delta_eta"))
      x@parameters$est[scaling] <- 1
      # x@parameters$std[scaling] <- NA
      x@parameters$par[scaling] <- 0
      # x@parameters$se[scaling] <- NA
      # x@parameters$p[scaling] <- NA
      # x@parameters$mi[scaling] <- NA
      # x@parameters$pmi[scaling] <- NA
      # x@parameters$mi_equal[scaling] <- NA
      # x@parameters$pmi_equal[scaling] <- NA
      x@parameters$fixed[scaling] <- TRUE
      x@parameters$identified[scaling] <- TRUE
      
      # Clear
      x@parameters <- clearpars(x@parameters, scaling)
    } else {
      # Set all first factor loadings equal to 1:
      for (i in unique(x@parameters$col[x@parameters$matrix == "lambda"])){
        firstLoading <- which(!x@parameters$fixed & x@parameters$matrix == "lambda" & x@parameters$col == i)[1]  
        x@parameters$est[firstLoading] <- 1
        x@parameters$par[firstLoading] <- 0
        x@parameters$fixed[firstLoading] <- TRUE
        x@parameters$identified[firstLoading] <- TRUE
        
        # Clear
        x@parameters <- clearpars(x@parameters, firstLoading)
      }
    }
    
    
    # Fix labels:
    x@parameters <- parRelabel(x@parameters)
    
  } else {
    # Number of equality constrains:
    cons <- x@parameters %>% group_by_("matrix","row","col") %>% summarize_(eq = ~!(all(fixed))&allTheSame(par))
    consPerMat <- cons %>% group_by_("matrix") %>% summarize_(n = ~sum(eq))
    
    nLat <- max(cons$col[cons$matrix == "lambda"])
    nMan <- max(cons$row[cons$matrix == "lambda"])
    
    ### LATENT MEANS ###
    # at least n_eta intercepts nead to be equal
    if (consPerMat$n[consPerMat$matrix == "tau"] >= nLat){
      means <- which(x@parameters$matrix %in% c("mu_eta") & x@parameters$group_id == 1)
      free <-  which(x@parameters$matrix %in% c("mu_eta") & x@parameters$group_id > 1)
    } else {
      means <- which(x@parameters$matrix %in% c("mu_eta"))
      free <- numeric(0)
    }
    
    # Constrain means:
    x@parameters$est[means] <- 0
    # x@parameters$std[means] <- NA
    x@parameters$par[means] <- 0
    # x@parameters$se[means] <- NA
    # x@parameters$p[means] <- NA
    # x@parameters$mi[means] <- NA
    # x@parameters$pmi[means] <- NA
    # x@parameters$mi_equal[means] <- NA
    # x@parameters$pmi_equal[means] <- NA
    x@parameters$fixed[means] <- TRUE
    x@parameters$identified[means] <- TRUE
    
    # Clear
    
    x@parameters <- clearpars(x@parameters, means)
    
    if (length(free) > 0){
      # x@parameters$std[free] <- NA
      x@parameters$par[free] <- max(x@parameters$par) + seq_along(free)
      # x@parameters$se[free] <- NA
      # x@parameters$p[free] <- NA
      # x@parameters$mi[free] <- NA
      # x@parameters$pmi[free] <- NA
      # x@parameters$mi_equal[free] <- NA
      # x@parameters$pmi_equal[free] <- NA
      x@parameters$fixed[free] <- FALSE
      x@parameters$identified[free] <- FALSE
      
      # Clear
      x@parameters <- clearpars(x@parameters, free)
    }
    
    
    if (type == "scaling"){
      ### SCALING ###
      # At least n_eta factor loadings need to be equal (FIXME: not sure about this...)
      if (consPerMat$n[consPerMat$matrix == "lambda"] >= nLat){
        scaling <- which(x@parameters$matrix %in% c("delta_eta") & x@parameters$group_id == 1)
        free <- which(x@parameters$matrix %in% c("delta_eta") & x@parameters$group_id > 1)
      } else {
        scaling <- which(x@parameters$matrix %in% c("delta_eta"))
        free <- numeric(0)
      }
      
      # Constrain scaling:
      # Set al scaling to 1:
      x@parameters$est[scaling] <- 1
      x@parameters$par[scaling] <- 0
      x@parameters$fixed[scaling] <- TRUE
      x@parameters$identified[scaling] <- TRUE
      x@parameters <- clearpars(x@parameters, free)
      
      if (length(free) > 0){
        x@parameters$par[free] <- max(x@parameters$par) + seq_along(free)
        x@parameters$fixed[free] <- FALSE
        x@parameters$identified[free] <- FALSE
        x@parameters <- clearpars(x@parameters, free)
      }
      
    } else {
      for (g in seq_len(nrow(x@sample@groups))){
        # Set all first factor loadings equal to 1:
        for (i in unique(x@parameters$col[x@parameters$matrix == "lambda"])){
  
          firstLoading <- which((!x@parameters$fixed | (x@parameters$fixed & x@parameters$est != 0)) & x@parameters$matrix == "lambda" & x@parameters$col == i & x@parameters$group_id == g)[1] 
          x@parameters$est[firstLoading] <- 1
          x@parameters$par[firstLoading] <- 0
          x@parameters$fixed[firstLoading] <- TRUE
          x@parameters$identified[firstLoading] <- TRUE
          
          # Clear
          x@parameters <- clearpars(x@parameters, firstLoading)
        }        
      }

    }   
    
  }
  

  
  # Fix labels:
  x@parameters <- parRelabel(x@parameters)
  # Return model:
  return(x)
}