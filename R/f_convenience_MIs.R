# psychonetrics parameter extraction:
MIs <- function(x,all = FALSE,matrices, sortby = c("mi","mi_equal"), top = 10,verbose = TRUE){
  sortby <- match.arg(sortby)
  if (missing(matrices)) matrices <- x@matrices$name
  
  # AWESOME HEADER!!!
  psychonetrics_print_logo()
  # cat(
  #   paste0("\t########################################\n",
  #          "\t## psychonetrics modification indices ##\n",
  #          "\t########################################\n\n"))
  # Obtain the parameter table:
  parTable <- x@parameters

  # filter only non-zero parameters and select only relevant columns:
  # parTable <- parTable %>% filter_(~fixed,~!is.na(mi)) %>% 
  #   select_("var1","op","var2","mi","pmi","mi_equal","pmi_equal","matrix","row","col","group")
  parTable <- parTable %>% filter_(~mi!=0) %>% 
    select_("var1","op","var2","est","mi","pmi","mi_equal","pmi_equal","matrix","row","col","group")
  
  # If nothing, return:
  if (nrow(parTable)==0){
    if (verbose) message("No modification indices to be reported...")
    return(NULL)
  }
  
  if (!all){
    # Display the top x:
    topx <- parTable[order(parTable[[sortby]],decreasing = TRUE),] %>% head(top)
    topx$mi <- goodNum(topx$mi)
    topx$pmi <- goodNum(topx$pmi)
    topx$mi_equal <- goodNum(topx$mi_equal)
    topx$pmi_equal <- goodNum(topx$pmi_equal)
    cat("Top ",top,"modification indices ",paste0("(ordered by ",sortby,")\n\n"))
    print.data.frame(topx, row.names=FALSE)
  } else {
    parTable <- parTable[order(parTable[[sortby]],decreasing = TRUE),] 
    # Make the entire table nice:
    parTable$mi <- goodNum(parTable$mi)
    parTable$pmi <- goodNum(parTable$pmi)
    parTable$mi_equal <- goodNum(parTable$mi_equal)
    parTable$pmi_equal <- goodNum(parTable$pmi_equal)
    
    
    # For each group:
    for (g in x@sample@groups$label){
      cat("\n Modification indices for group",g,"\n")
      # for each matrix:
      for (mat in unique(parTable$matrix[parTable$group == g])){
        if (x@matrices$diagonal[x@matrices$name==mat]){
          brackets <- "(diagonal)"
        } else  if (x@matrices$symmetrical[x@matrices$name==mat]){
          brackets <- "(symmetric)"
        }  else {
          brackets <- "" 
        }
        subTable <- parTable %>% filter_(~group == g,~matrix==mat) %>% select_(~-matrix,~-group) 
        # subTable <- subTable[order(subTable[[sortby]],decreasing = TRUE),] 
        if (nrow(subTable) > 0){
          cat("\n\t- ",mat,brackets,"\n")
          print.data.frame(subTable, row.names=FALSE)          
        }
      }
    }
    
  }
  cat("\nNote: mi_equal = modification index if parameter is added in all groups (constrained equal)")
  
  
  invisible(x@parameters %>% filter_(~fixed,~!is.na(mi)) %>% 
              select_("var1","op","var2","mi","pmi","mi_equal","pmi_equal","matrix","row","col","group"))
}
