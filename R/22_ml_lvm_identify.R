
# lvm identifier:
identify_ml_lvm <- function(x){
  # Idenitfy type:
  type <- x@identification

  
  # Always set beta diagonal to zero:
  
  betaDiag <- which(grepl("beta",x@parameters$matrix) & x@parameters$row == x@parameters$col)
  x@parameters$est[betaDiag] <- 0
  x@parameters$par[betaDiag] <- 0
  x@parameters$fixed[betaDiag] <- TRUE
  x@parameters$identified[betaDiag] <- TRUE
  x@parameters <- clearpars(x@parameters, betaDiag)
  
  # Always set residual (co)variance of single indicator items to zero:
  nIndicators <- x@parameters %>% filter(matrix == "lambda") %>% group_by(.data[['var2_id']]) %>% 
    summarize(nLats = length(unique(.data[['var1_id']][!.data[['fixed']]])))
  
  if (any(nIndicators$nLats == 1)){
    for (inds in unique(x@parameters$var1_id[x@parameters$var2_id%in%nIndicators$var2_id[nIndicators$nLats == 1] & x@parameters$matrix == "lambda"])){
      # Which to constrain?
      cons <- x@parameters$var1_id == inds & 
        x@parameters$matrix %in% c("sigma_epsilon_within","omega_epsilon_within","delta_epsilon_within","lowertri_epsilon_within","kappa_epsilon_within",
                                   "sigma_epsilon_between","omega_epsilon_between","delta_epsilon_between","lowertri_epsilon_between","kappa_epsilon_between")
      
      x@parameters$est[cons] <- 0
      x@parameters$par[cons] <- 0
      x@parameters$fixed[cons] <- TRUE
      x@parameters$identified[cons] <- TRUE
      x@parameters <- clearpars(x@parameters, cons)
    }
  }
  
  
  # Single group is easy:
  if (nrow(x@sample@groups) == 1){
    
    # Set all latent intercepts to zero:
    means <- which(x@parameters$matrix %in% c("mu_eta"))
    x@parameters$est[means] <- 0
    x@parameters$par[means] <- 0
    x@parameters$fixed[means] <- TRUE
    x@parameters$identified[means] <- TRUE
    
    x@parameters <- clearpars(x@parameters, means)
    
    
    # variance ifentification:
    if (type == "variance"){
      # 
      # Within-subject model
      mat <- switch(
       x@types$within_latent,
       "cov" = "sigma_zeta_within",
       "prec" = "kappa_zeta_within",
       "ggm" = "delta_zeta_within",
       "chol" = "lowertri_zeta_within"
      )

      # Set all latent variances to 1:
      vars <- which(x@parameters$matrix == mat & x@parameters$row == x@parameters$col)
      x@parameters$est[vars] <- 1
      x@parameters$par[vars] <- 0
      x@parameters$fixed[vars] <- TRUE
      x@parameters$identified[vars] <- TRUE

      # Clear
      x@parameters <- clearpars(x@parameters, vars)
      # 
      # Between-subject model
      # mat <- switch(
      #   x@types$between_latent,
      #   "cov" = "sigma_zeta_between",
      #   "prec" = "kappa_zeta_between",
      #   "ggm" = "delta_zeta_between",
      #   "chol" = "lowertri_zeta_between"
      # )
      # 
      # # Set all latent variances to 1:
      # vars <- which(x@parameters$matrix == mat & x@parameters$row == x@parameters$col)
      # x@parameters$est[vars] <- 1
      # x@parameters$par[vars] <- 0
      # x@parameters$fixed[vars] <- TRUE
      # x@parameters$identified[vars] <- TRUE
      # 
      # # Clear
      # x@parameters <- clearpars(x@parameters, vars)
      
    } else {
      # Within-subject part:
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
      # 
      # # Between-subject part:
      # # Set all first factor loadings equal to 1:
      # for (i in unique(x@parameters$col[x@parameters$matrix == "lambda"])){
      #   # firstLoading <- which(!x@parameters$fixed & x@parameters$matrix == "lambda" & x@parameters$col == i)[1]  
      #   firstLoading <- which((!x@parameters$fixed | (x@parameters$fixed & x@parameters$est != 0)) & x@parameters$matrix == "lambda" & x@parameters$col == i)[1] 
      #   
      #   x@parameters$est[firstLoading] <- 1
      #   x@parameters$par[firstLoading] <- 0
      #   x@parameters$fixed[firstLoading] <- TRUE
      #   x@parameters$identified[firstLoading] <- TRUE
      #   
      #   # Clear
      #   x@parameters <- clearpars(x@parameters, firstLoading)
      # }
    }
    
    
    # Fix labels:
    x@parameters <- parRelabel(x@parameters)
    
  } else {
    # Number of equality constrains:
    cons <- x@parameters %>% group_by(.data[["matrix"]],.data[["row"]],.data[["col"]]) %>% summarize(eq = !(all(.data[['fixed']]))&allTheSame(.data[['par']]))
    consPerMat <- cons %>% group_by(.data[["matrix"]]) %>% summarize(n = sum(.data[['eq']]))
    
    nLat <- max(cons$col[cons$matrix == "lambda"])
    # nLat <- max(cons$col[cons$matrix == "lambda"])
    nMan <- max(cons$row[cons$matrix == "lambda"])
    
    ### LATENT MEANS ###
    # at least n_eta intercepts nead to be equal
    if (consPerMat$n[consPerMat$matrix == "nu"] >= nLat){
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
    
    
    if (type == "variance"){
      

      ### variance ###
      # At least n_eta factor loadings need to be equal (FIXME: not sure about this...)
      # Between
      # mat <- switch(
      #   x@types$between_latent,
      #   "cov" = "sigma_zeta_between",
      #   "prec" = "kappa_zeta_between",
      #   "ggm" = "delta_zeta_between",
      #   "chol" = "lowertri_zeta_between"
      # )
      # 
      # 
      # if (consPerMat$n[consPerMat$matrix == "lambda"] >= nLat){
      #   variance <- which(x@parameters$matrix == mat & x@parameters$group_id == 1 & x@parameters$row == x@parameters$col)
      #   free <- which(x@parameters$matrix == mat & x@parameters$group_id > 1 & x@parameters$row == x@parameters$col)
      # } else {
      #   variance <- which(x@parameters$matrix == mat &  x@parameters$row == x@parameters$col)
      #   free <- numeric(0)
      # }
      # 
      # # Constrain variance:
      # # Set al variance to 1:
      # x@parameters$est[variance] <- 1
      # x@parameters$par[variance] <- 0
      # x@parameters$fixed[variance] <- TRUE
      # x@parameters$identified[variance] <- TRUE
      # x@parameters <- clearpars(x@parameters, free)
      # 
      # if (length(free) > 0){
      #   x@parameters$par[free] <- max(x@parameters$par) + seq_along(free)
      #   x@parameters$fixed[free] <- FALSE
      #   x@parameters$identified[free] <- FALSE
      #   x@parameters <- clearpars(x@parameters, free)
      # }
      
      # Within
      mat <- switch(
        x@types$within_latent,
        "cov" = "sigma_zeta_within",
        "prec" = "kappa_zeta_within",
        "ggm" = "delta_zeta_within",
        "chol" = "lowertri_zeta_within"
      )
      
      if (consPerMat$n[consPerMat$matrix == "lambda"] >= nLat){
        variance <- which(x@parameters$matrix == mat & x@parameters$group_id == 1 & x@parameters$row == x@parameters$col)
        free <- which(x@parameters$matrix == mat & x@parameters$group_id > 1 & x@parameters$row == x@parameters$col)
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
        # # Set all first factor loadings equal to 1:
        # for (i in unique(x@parameters$col[x@parameters$matrix == "lambda"])){
        #   
        #   firstLoading <- which((!x@parameters$fixed | (x@parameters$fixed & x@parameters$est != 0)) & x@parameters$matrix == "lambda" & x@parameters$col == i & x@parameters$group_id == g)[1] 
        #   x@parameters$est[firstLoading] <- 1
        #   x@parameters$par[firstLoading] <- 0
        #   x@parameters$fixed[firstLoading] <- TRUE
        #   x@parameters$identified[firstLoading] <- TRUE
        #   
        #   # Clear
        #   x@parameters <- clearpars(x@parameters, firstLoading)
        # }
        # 
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