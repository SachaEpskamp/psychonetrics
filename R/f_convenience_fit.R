# psychonetrics parameter extraction:
fit <- function(x){
  
  if (is(x,"psychonetrics_bootstrap")){
    fittable <- x@fitmeasures[,c("measure","min","q2.5","q97.5","max")]
    for (j in 2:ncol(fittable)){
      fittable[[j]] <- goodNum2(fittable[[j]])
    }
    
    print.data.frame(fittable,row.names=FALSE)
    invisible(x@fitmeasures)
  } else {
    
    # Bootstrap warning:
    if (x@sample@bootstrap){
      boot_warning()
    }
    
    
    if (!x@computed){
      stop("Model has not yet been computed. Use 'runmodel'.")
    }
    if (is.null(x@fitmeasures)){
      stop("No fit measures computed yet. Use 'addfit'.")
    }
    
    # AWESOME HEADER!!!
    # cat(
    #   paste0("\t################################\n",
    #          "\t## psychonetrics fit measures ##\n",
    #          "\t################################\n\n"))
    # psychonetrics_print_logo() 
    # No awesome header :(
    
    df <- data.frame(
      Measure = names(x@fitmeasures),
      Value = goodNum2(unlist(x@fitmeasures))
    )
    print.data.frame(df,row.names=FALSE)
    
    # Numeric DF to return:
    
    dfnum <- data.frame(
      Measure = names(x@fitmeasures),
      Value = unlist(x@fitmeasures)
    )
    invisible(dfnum)
  }
}
