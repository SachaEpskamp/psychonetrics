
# lvm identifier:
identify_lvm <- function(x){
  # Idenitfy type:
  type <- x@identification

  
  # Always set beta diagonal to zero:
  betaDiag <- which(x@parameters$matrix == "beta" & x@parameters$row == x@parameters$col)
  x@parameters$est[betaDiag] <- 0
  x@parameters$par[betaDiag] <- 0
  x@parameters$fixed[betaDiag] <- TRUE
  x@parameters$identified[betaDiag] <- TRUE
  x@parameters <- clearpars(x@parameters, betaDiag)
  
  # Always set residual (co)variance of single indicator items to zero:
  nIndicators <- x@parameters %>% filter(.data[['matrix']] == "lambda") %>% group_by(.data[['var2_id']]) %>% 
    summarize(nLats = length(unique(.data[['var1_id']][!.data[['fixed']] | .data[['est']] != 0])))

  if (any(nIndicators$nLats == 1)){
    for (inds in unique(x@parameters$var1_id[x@parameters$var2_id%in%nIndicators$var2_id[nIndicators$nLats == 1] & x@parameters$matrix == "lambda" & !x@parameters$fixed])){
      # Which to constrain?
      cons <- x@parameters$var1_id == inds & 
        x@parameters$matrix %in% c("sigma_epsilon","omega_epsilon","delta_epsilon","lowertri_epsilon","kappa_epsilon")
      
      x@parameters$est[cons] <- 0
      x@parameters$par[cons] <- 0
      x@parameters$fixed[cons] <- TRUE
      x@parameters$identified[cons] <- TRUE
      x@parameters <- clearpars(x@parameters, cons)
    }
  }
  
  # Single group is easy:
  if (nrow(x@sample@groups) == 1){
    
    # How many values of nu are fixed?
    nConsInTau <- sum(x@parameters$fixed[x@parameters$matrix == "nu"])
    
    # Set all latent intercepts to zero if there are not enough constrains in nu:
    means <- which(x@parameters$matrix %in% c("nu_eta"))
    if (nConsInTau <= length(means)){
      
      x@parameters$est[means] <- 0
      x@parameters$par[means] <- 0
      x@parameters$fixed[means] <- TRUE
      x@parameters$identified[means] <- TRUE
      
      x@parameters <- clearpars(x@parameters, means)
    }

    
    
    # variance ifentification:
    if (type == "variance"){

      mat <- switch(
       x@types$latent,
       "cov" = "sigma_zeta",
       "prec" = "kappa_zeta",
       "ggm" = "delta_zeta",
       "chol" = "lowertri_zeta"
      )
      
      # Set all latent variances to 1:
      vars <- which(x@parameters$matrix == mat & x@parameters$row == x@parameters$col)
      x@parameters$est[vars] <- 1
      x@parameters$par[vars] <- 0
      x@parameters$fixed[vars] <- TRUE
      x@parameters$identified[vars] <- TRUE
      
      # Clear
      x@parameters <- clearpars(x@parameters, vars)
    } else {
      
      # Set all first factor loadings equal to 1:
      for (i in unique(x@parameters$col[x@parameters$matrix == "lambda"])){
        # firstLoading <- which(!x@parameters$fixed & x@parameters$matrix == "lambda" & x@parameters$col == i)[1]  
        firstLoading <- which((!x@parameters$fixed | (x@parameters$fixed & x@parameters$est != 0)) & x@parameters$matrix == "lambda" & x@parameters$col == i)[1] 
        
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
    cons <- x@parameters %>% group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% summarize(
      eq = !(all(.data[['fixed']]))&allTheSame(.data[['par']]))
    consPerMat <- cons %>% group_by(.data[["matrix"]]) %>% summarize(n = sum(.data[['eq']]))
    
    nLat <- max(cons$col[cons$matrix == "lambda"])
    nMan <- max(cons$row[cons$matrix == "lambda"])
    
    # FIXME:
    nFixed <- NULL
    fixedingroup1 <- x@parameters %>%  group_by(.data[["matrix"]],.data[["group"]]) %>% summarize(
     nFixed = sum(.data[['fixed']])) %>% group_by(.data[["matrix"]]) %>% summarize(nFixed = max(.data[['nFixed']]))
     
    nNuFix <- fixedingroup1$nFixed[fixedingroup1$matrix == "nu"]
    
    ### LATENT MEANS ###
    # at least n_eta intercepts nead to be equal
    if (!is.null(nNuFix) && nNuFix == 0){
      if (consPerMat$n[consPerMat$matrix == "nu"] >= nLat){
        means <- which(x@parameters$matrix %in% c("nu_eta") & x@parameters$group_id == 1)
        free <-  which(x@parameters$matrix %in% c("nu_eta") & x@parameters$group_id > 1 & !(x@parameters$fixed & !x@parameters$identified))
      } else {
        means <- which(x@parameters$matrix %in% c("nu_eta"))
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
      
    }

    
    if (type == "variance"){
      
      mat <- switch(
        x@types$latent,
        "cov" = "sigma_zeta",
        "prec" = "kappa_zeta",
        "ggm" = "delta_zeta",
        "chol" = "lowertri_zeta"
      )
      
      ### variance ###
      # At least n_eta factor loadings need to be equal (FIXME: not sure about this...)
      if (consPerMat$n[consPerMat$matrix == "lambda"] >= nLat){
        variance <- which(x@parameters$matrix == mat & x@parameters$group_id == 1 & x@parameters$row == x@parameters$col)
        free <- which(x@parameters$matrix == mat & x@parameters$group_id > 1 & x@parameters$row == x@parameters$col & x@parameters$identified)
        # free <- which(x@parameters$matrix == mat & x@parameters$group_id > 1 & x@parameters$row == x@parameters$col)
      } else {
        variance <- which(x@parameters$matrix == mat &  x@parameters$row == x@parameters$col)
        free <- numeric(0)
      }
      
      # Constrain variance:
      # Set al variance to 1:
      x@parameters$est[variance] <- 1
      x@parameters$par[variance] <- 0
      x@parameters$fixed[variance] <- TRUE
      x@parameters$identified[variance] <- TRUE
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