# psychonetrics parameter extraction:
parameters <- function(x){
  
  # FIXME: if the model is aggregate bootstraps run different function instead:
  if (is(x,"psychonetrics_bootstrap")){
    return(parameters_bootstrap(x))
  }
  
  # Error if not psychonetrics:
  stopifnot(is(x,"psychonetrics"))
  
  
  # Bootstrap warning:
  if (x@sample@bootstrap){
    boot_warning()
  }
  

  # # AWESOME HEADER!!!
  # cat(
  #   paste0("\t##################################\n",
  #          "\t## psychonetrics parameter list ##\n",
  #          "\t##################################\n\n"))
  # psychonetrics_print_logo()
  # No awesome header :(
  # Obtain the parameter table:
  parTable <- x@parameters
  
  # filter only non-zero parameters and select only relevant columns:
  # FIXME: std not yet implemented, so remove now:
  
  if (x@estimator == "PML") {
    cols <- c("var1","op","var2","est","penalty_lambda","matrix","row","col","group","par")
    boots <- FALSE
  } else if (!all(is.na(x@parameters$se_boot))){
    cols <- c("var1","op","var2","est","se","p","se_boot","p_boot","matrix","row","col","group","par")
    boots <- TRUE
  } else {
    cols <- c("var1","op","var2","est","se","p","matrix","row","col","group","par")
    boots <- FALSE
  }
  parTable <- parTable %>% filter(drop(!.data[['fixed']]|.data[['est']] != 0)) %>%
    select(all_of(cols))

  # Make var2 nicer:
  parTable$var2 <- ifelse(is.na(parTable$var2),"",parTable$var2)

  # Pretty numbers:
  parTable$est <- goodNum2(parTable$est)
  if (x@estimator == "PML") {
    parTable$penalty_lambda <- goodNum2(parTable$penalty_lambda)
  } else {
    parTable$se <- goodNum(parTable$se)
    parTable$p <- goodNum(parTable$p)
  }

  if (boots){
    parTable$se_boot <- goodNum(parTable$se_boot)
    parTable$p_boot <- goodNum(parTable$p_boot)
  }
  
  
  # if not computed, remove est, se and p:
  if (!x@computed){
    warning("Model has not been computed! Not returning estimates, standard errors and p-values. Use runmodel() to compute the model.")
    parTable <- parTable %>% select(-.data[['est']],-.data[['se']],-.data[['p']])
  }
  
  # For each group:
  for (g in x@sample@groups$label){
    cat("\n Parameters for group",g)
    # for each matrix:
    for (mat in unique(parTable$matrix[parTable$group == g])){
      if (x@matrices$diagonal[x@matrices$name==mat]){
        brackets <- "(diagonal)"
      } else  if (x@matrices$symmetrical[x@matrices$name==mat]){
        brackets <- "(symmetric)"
      }  else {
        brackets <- "" 
      }
      cat("\n\t- ",mat,brackets,"\n")
      print.data.frame(parTable %>% filter(.data[['group']] == g,.data[['matrix']]==mat) %>% select(-.data[['matrix']],-.data[['group']]), row.names=FALSE)
    }
  }
  
  invisible(x@parameters)
}

# Bootstrap version:
parameters_bootstrap <- function(x){
  
  parTable <- x@parameters
  
  # filter only non-zero parameters and select only relevant columns:
  # FIXME: std not yet implemented, so remove now:
  # 
  # if (!all(is.na(x@parameters$se_boot))){
  #   cols <- c("var1","op","var2","est","se","p","se_boot","p_boot","matrix","row","col","group","par")
  #   boots <- TRUE
  # } else {
  #   cols <- c("var1","op","var2","est","se","p","matrix","row","col","group","par")
  #   boots <- FALSE
  # }
  
  # columns for all pars:
  cols_incl0 <- c("var1", "op", "var2", "est_sample", "se_boot", 
             "q2.5", "q97.5", "q99", 
            "matrix", "row", 
            "col", "group")

  # columns for non-zero pars:
  cols_non0 <- c("var1", "op", "var2", "est_sample", "se_boot_non0", "prop_non0", "prop_non0_pos", "prop_non0_neg", 
             "q2.5_non0", "q97.5_non0", 
            "matrix", "row", 
            "col", "group")
  
  
  # Make var2 nicer:
  parTable$var2 <- ifelse(is.na(parTable$var2),"",parTable$var2)
  
  # Pretty numbers:
  if (!any(is.na(parTable$fixed_sample))){
    parTable$est_sample <- ifelse(parTable$fixed_sample,"",goodNum2(parTable$est_sample))
  } else {
    parTable$est_sample <- goodNum2(parTable$est_sample)
  }
  
  parTable$avg <- goodNum2(parTable$avg)
  parTable$min <- goodNum2(parTable$min)
  parTable$q1 <- goodNum2(parTable$q1)
  parTable$q2.5 <- goodNum2(parTable$q2.5)
  parTable$q5 <- goodNum2(parTable$q5)
  parTable$median <- goodNum2(parTable$median)
  parTable$q95 <- goodNum2(parTable$q95)
  parTable$q97.5 <- goodNum2(parTable$q97.5)
  parTable$q99 <- goodNum2(parTable$q99)
  parTable$max <- goodNum2(parTable$max)
  
  parTable$avg_non0 <- goodNum2(parTable$avg_non0)
  parTable$min_non0 <- goodNum2(parTable$min_non0)
  parTable$q1_non0 <- goodNum2(parTable$q1_non0)
  parTable$q2.5_non0 <- goodNum2(parTable$q2.5_non0)
  parTable$q5_non0 <- goodNum2(parTable$q5_non0)
  parTable$median_non0 <- goodNum2(parTable$median_non0)
  parTable$q95_non0 <- goodNum2(parTable$q95_non0)
  parTable$q97.5_non0 <- goodNum2(parTable$q97.5_non0)
  parTable$q99_non0 <- goodNum2(parTable$q99_non0)
  parTable$max_non0 <- goodNum2(parTable$max_non0)
  
  parTable$se_boot <- goodNum(parTable$se_boot)
  parTable$se_boot_non0 <- goodNum(parTable$se_boot_non0)
  parTable$prop_non0 <- goodNum3(parTable$prop_non0)
  parTable$prop_non0_pos <- goodNum3(parTable$prop_non0_pos)
  parTable$prop_non0_neg <- goodNum3(parTable$prop_non0_neg)
  

  # For each group:
  for (g in unique(parTable$group)){
    cat("\n Parameters for group",g)
    # for each matrix:
    for (mat in unique(parTable$matrix[parTable$group == g])){
      if (x@matrices$diagonal[x@matrices$name==mat]){
        brackets <- "(diagonal)"
      } else  if (x@matrices$symmetrical[x@matrices$name==mat]){
        brackets <- "(symmetric)"
      }  else {
        brackets <- "" 
      }
      cat("\n\t- ",mat,brackets,"-- including times parameter was zero\n")
      print.data.frame(parTable %>% select(all_of(cols_incl0)) %>% filter(.data[['group']] == g,.data[['matrix']]==mat) %>% select(-.data[['matrix']],-.data[['group']]), row.names=FALSE)
      
      cat("\n\t- ",mat,brackets,"-- *not* including times parameter was zero\n")
      df_print <- parTable %>% select(all_of(cols_non0)) %>% filter(.data[['group']] == g,.data[['matrix']]==mat) %>% select(-.data[['matrix']],-.data[['group']])
      names(df_print)[names(df_print)!="prop_non0"] <- gsub("_non0","",names(df_print)[names(df_print)!="prop_non0"])
      print.data.frame(df_print, row.names=FALSE)
    }
  }
  
  invisible(x@parameters)
}
