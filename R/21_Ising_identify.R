
# Ising identifier:
identify_Ising <- function(x){

  # Which parameters are the beta parameters?
  betas <- which(x@parameters$matrix == "beta")
  
  # In a single group analysis, beta is always fixed and identified:
  if (nrow(x@sample@groups) == 1){
    x@parameters$par[betas] <- 0
    x@parameters$fixed[betas] <- TRUE
    x@parameters$identified[betas] <- TRUE
  } else {
    # If there are equality constrains across groups in omega or tau, free all temperature
    
    # Number of equality constrains per matrix:
    cons <- x@parameters %>% group_by_("matrix","row","col") %>% summarize_(eq = ~!(all(fixed))&allTheSame(par))
    consPerMat <- cons %>% group_by_("matrix") %>% summarize_(n = ~sum(eq))
    
    # at least 1 intercepts nead to be equal
    if (sum(consPerMat$n[consPerMat$matrix %in% c("omega","tau")]) >= 1){
      fix <- which(x@parameters$matrix %in% c("beta") & x@parameters$group_id == 1)
      free <-  which(x@parameters$matrix %in% c("beta") & x@parameters$group_id > 1 & !(x@parameters$fixed & !x@parameters$identified))
    } else {
      fix <- which(x@parameters$matrix %in% c("beta"))
      free <- numeric(0)
    }
    # Constrain means:
    x@parameters$est[fix] <- 1
    # x@parameters$std[means] <- NA
    x@parameters$par[fix] <- 0
    # x@parameters$se[means] <- NA
    # x@parameters$p[means] <- NA
    # x@parameters$mi[means] <- NA
    # x@parameters$pmi[means] <- NA
    # x@parameters$mi_equal[means] <- NA
    # x@parameters$pmi_equal[means] <- NA
    x@parameters$fixed[fix] <- TRUE
    x@parameters$identified[fix] <- TRUE
    
    # Clear
    
    x@parameters <- clearpars(x@parameters, fix)
    
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
    
    
    
  }
  
  # Fix labels:
  x@parameters <- parRelabel(x@parameters)
  

  
  # Return model:
  return(x)
}