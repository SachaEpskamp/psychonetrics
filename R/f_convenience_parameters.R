# psychonetrics parameter extraction:
parameters <- function(x){
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
  
  if (!all(is.na(x@parameters$se_boot))){
    cols <- c("var1","op","var2","est","se","p","se_boot","p_boot","matrix","row","col","group","par")
    boots <- TRUE
  } else {
    cols <- c("var1","op","var2","est","se","p","matrix","row","col","group","par")
    boots <- FALSE
  }
  parTable <- parTable %>% filter(!.data[['fixed']]|.data[['est']] != 0) %>% 
    select(all_of(cols))
  
  # Make var2 nicer:
  parTable$var2 <- ifelse(is.na(parTable$var2),"",parTable$var2)
  
  # Pretty numbers:
  parTable$est <- goodNum2(parTable$est)
  parTable$se <- goodNum(parTable$se)
  parTable$p <- goodNum(parTable$p)
  
  if (boots){
    parTable$se_boot <- goodNum(parTable$se_boot)
    parTable$p_boot <- goodNum(parTable$p_boot)
  }
  
  
  # if not computed, remove est, se and p:
  if (!x@computed){
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
