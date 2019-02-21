# psychonetrics parameter extraction:
parameters <- function(x){
  # # AWESOME HEADER!!!
  # cat(
  #   paste0("\t##################################\n",
  #          "\t## psychonetrics parameter list ##\n",
  #          "\t##################################\n\n"))
  psychonetrics_print_logo()
  # Obtain the parameter table:
  parTable <- x@parameters
  
  # filter only non-zero parameters and select only relevant columns:
  # FIXME: std not yet implemented, so remove now:
  parTable <- parTable %>% filter_(~!fixed|est != 0) %>% 
    select_("var1","op","var2","est","se","p","matrix","row","col","group","par")
  
  # Make var2 nicer:
  parTable$var2 <- ifelse(is.na(parTable$var2),"",parTable$var2)
  
  # Pretty numbers:
  parTable$est <- goodNum2(parTable$est)
  parTable$se <- goodNum(parTable$se)
  parTable$p <- goodNum(parTable$p)
  
  # if not computed, remove est, se and p:
  if (!x@computed){
    parTable <- parTable %>% select_(~-est,~-se,~-p)
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
      print.data.frame(parTable %>% filter_(~group == g,~matrix==mat) %>% select_(~-matrix,~-group), row.names=FALSE)
    }
  }
  

}