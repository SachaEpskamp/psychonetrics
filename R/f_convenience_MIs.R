MIs <- function(x, all = FALSE, matrices, type = c("normal","equal","free"), top = 10,verbose = TRUE, nonZero = FALSE){
  # Check if equal and free are needed:
  if (any(type != "normal")){
    if (nrow(x@sample@groups) == 1){
      type <- "normal"
    }
  }
  
  Results <- list()
  # Print the tables:
  for (t in seq_along(type)){
    Results[[t]] <- MIs_inner(x, all=all, matrices = matrices, type = type[t], top = top, verbose = verbose, nonZero = nonZero)
  }
  
  if (length(Results)==1){
    invisible(Results[[1]])
  } else {
    invisible(Results)
  }
}

# psychonetrics parameter extraction:
MIs_inner <- function(x,all = FALSE, matrices, type = c("normal","equal","free"), top = 10,verbose = TRUE, nonZero = FALSE){
  # sortby <- match.arg(sortby)
  if (missing(matrices)) matrices <- x@matrices$name
  
  # AWESOME HEADER!!!
  # psychonetrics_print_logo()
  # No awesome header :(
  # cat(
  #   paste0("\t########################################\n",
  #          "\t## psychonetrics modification indices ##\n",
  #          "\t########################################\n\n"))
  # Obtain the parameter table:
  parTable <- x@parameters
  
  # Which columns to use?
  micol <- switch(type,
                  "normal" = "mi",
                  "equal" = "mi_equal",
                  "free" = "mi_free"
                  )
  pcol <- switch(type,
                  "normal" = "pmi",
                  "equal" = "pmi_equal",
                  "free" = "pmi_free"
  )
  epccol <- switch(type,
                 "normal" = "epc",
                 "equal" = "epc_equal",
                 "free" = "epc_free"
  )

  # filter only non-zero parameters and select only relevant columns:
  parTable <- parTable %>%  filter(.data[['matrix']] %in% matrices) %>% 
    select(.data[["var1"]],.data[["op"]],.data[["var2"]],.data[["est"]],.data[[micol]],.data[[pcol]],.data[[epccol]],.data[["matrix"]],.data[["row"]],.data[["col"]],.data[["group"]],.data[["group_id"]])
  
  # nonZero:
  if (nonZero){
    parTable <- parTable %>% filter(.data[['est']]!=0)
  }
  
  # If nothing, return:
  if (nrow(parTable)==0){
    if (verbose) message("No modification indices to be reported...")
    return(NULL)
  }
  
  if (!all){
    # Display the top x:
    topx <- parTable[order(parTable[[micol]],decreasing = TRUE),] %>% head(top)
    topx[[micol]] <- goodNum(topx[[micol]])
    topx[[pcol]] <- goodNum(topx[[pcol]])
    topx[[epccol]] <- goodNum2(topx[[epccol]])
   
    # cat("Top ",top,"modification indices ",paste0("(ordered by ",sortby,")\n\n"))
    message <- switch(type,
              "normal" = paste0("\nTop ",top," modification indices:\n\n"),
              "equal" = paste0("\nTop ",top," group-constrained modification indices:\n\n"),
              "free" = paste0("\nTop ",top," equality-free modification indices:\n\n"))
    cat(message)
    print.data.frame(topx, row.names=FALSE)
    
  } else {
    parTable <- parTable[order(parTable[[micol]],decreasing = TRUE),] 
    # Make the entire table nice:
    parTable[[micol]] <- goodNum(parTable[[micol]])
    parTable[[pcol]] <- goodNum(parTable[[pcol]])
    parTable[[epccol]] <- goodNum2(parTable[[epccol]])

    # For each group:
    
    
    for (g in x@sample@groups$id){
      message <- switch(type,
                        "normal" = paste0("\nModification indices for group:",x@sample@groups$label[g]),
                        "equal" = paste0("\nGroup-constrained modification indicesfor group:",x@sample@groups$label[g]),
                        "free" = paste0("\nEquality-free modification indices for group:",x@sample@groups$label[g]))
      cat(message)
      
      # cat(paste0("\n\nGroup ",x@sample@groups$label[g]))
      
      # for each matrix:
      # for (mat in unique(parTable$matrix[parTable$group_id == g])){
      for (mat in matrices){
        if (x@matrices$diagonal[x@matrices$name==mat]){
          brackets <- "(diagonal)"
        } else  if (x@matrices$symmetrical[x@matrices$name==mat]){
          brackets <- "(symmetric)"
        }  else {
          brackets <- "" 
        }
        subTable <- parTable %>% filter(.data[['group_id']] == g,.data[['matrix']]==mat) %>% select(-.data[['matrix']],-.data[['group']],-.data[['group_id']]) 
        # subTable <- subTable[order(subTable[[sortby]],decreasing = TRUE),] 
        if (nrow(subTable) > 0){
          cat("\n\t- ",mat,brackets,"\n")
          print.data.frame(subTable, row.names=FALSE)          
        }
      }
    }
    
  }
  

  return(parTable)
}
