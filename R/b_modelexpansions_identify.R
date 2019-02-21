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
  # Single group is easy:
  if (nrow(x@sample@groups) == 1){
    
    # Set all means to zero:
    means <- which(x@parameters$matrix %in% c("mu_eta"))
    x@parameters$est[means] <- 0
    x@parameters$std[means] <- NA
    x@parameters$par[means] <- 0
    x@parameters$se[means] <- NA
    x@parameters$p[means] <- NA
    x@parameters$mi[means] <- NA
    x@parameters$pmi[means] <- NA
    x@parameters$mi_equal[means] <- NA
    x@parameters$pmi_equal[means] <- NA
    x@parameters$fixed[means] <- TRUE
    x@parameters$identified[means] <- TRUE
    
    # Set al scaling to 1:
    scaling <- which(x@parameters$matrix %in% c("delta_eta"))
    x@parameters$est[scaling] <- 1
    x@parameters$std[scaling] <- NA
    x@parameters$par[scaling] <- 0
    x@parameters$se[scaling] <- NA
    x@parameters$p[scaling] <- NA
    x@parameters$mi[scaling] <- NA
    x@parameters$pmi[scaling] <- NA
    x@parameters$mi_equal[scaling] <- NA
    x@parameters$pmi_equal[scaling] <- NA
    x@parameters$fixed[scaling] <- TRUE
    x@parameters$identified[scaling] <- TRUE

    # Fix labels:
    x@parameters <- parRelabel(x@parameters)
    
  } else {
    stop("Multi-group LNM identification is not yet implemented. Please send a dissapointed mail to the maintainer!")
  }
  
  # Return model:
  return(x)
}