# Compare function for psychonetrics models:
compare <- function(...){
  # Obtain dots:
  dots <- list(...)
  
  # if anything is not a psychonetrics object, stop:
  classes <- sapply(dots, class)
  if (!all(classes == "psychonetrics")){
    stop("psychonetrics::compare(...) must only contain psychonetrics objects")
  }
  
  # If not all computed, stop:
  if (!all(sapply(dots,function(x)x@computed))){
    stop("Not all models are computed.")
  }
  
  # If names are null, add names:
  if (is.null(names(dots))){
    names(dots) <- paste0("Model ",seq_along(dots))
  }
  
  # If there are not at least two models, compare to saturated:
  if (length(dots)==1){
    stop("psychonetrics::compare(...) requires at least two models")  
    # dots <- c(list(saturated = dots[[1]]@baseline_saturated$saturated),dots,list(baseline = dots[[1]]@baseline_saturated$saturated))
  }
  
  # Obtain the fitmeasures of each model:
  fits <- lapply(dots, function(x) x@fitmeasures)

  # Create first table:
  Tab <- data.frame(
    model = names(dots),
    DF = sapply(fits,"[[","df"),
    AIC = sapply(fits,"[[","aic.ll"),
    BIC = sapply(fits,"[[","bic"),
    RMSEA = sapply(fits,"[[","rmsea"),
    Chisq = sapply(fits,"[[","chisq"),
    Chisq_diff = NA,
    DF_diff = NA,
    p_value = NA
  )
  
  # Arrange table by Df:
  Tab <- Tab %>% arrange(.data[['DF']])
  
  # Compute chis difference:
  Tab$Chisq_diff <- c(NA,abs(diff(Tab$Chisq)))
  
  # Compute DF diff:
  Tab$DF_diff <- c(NA,abs(diff(Tab$DF)))
  
  # Compute p:
  Tab$p_value <- pchisq(Tab$Chisq_diff,Tab$DF_diff,lower.tail=FALSE)
  
  # Set saturated chisquare to NA:
  if (any(Tab$model == "saturated"))
  Tab$Chisq[Tab$model == "saturated"] <- NA
  
  class(Tab) <- c("psychonetrics_compare","data.frame")
  
  return(Tab)
}

# Nice print function:
print.psychonetrics_compare <- function(x,...){

  x$AIC <- goodNum2(x$AIC)
  x$BIC <- goodNum2(x$BIC)
  x$RMSEA <- goodNum2(x$RMSEA)
  x$Chisq <- goodNum2(x$Chisq)
  x$Chisq_diff <- goodNum2(x$Chisq_diff)
  x$DF_diff <- goodNumInt(x$DF_diff)
  x$p_value <- goodNum(x$p_value)
  
  
  # Make all numbers nicer
  # for (i in 5:ncol(x)){
    # x[,i] <- goodNum(x[,i])
  # }
  
  # Output something:
  # cat(
  #   paste0("\t\t####################################\n",
  #          "\t\t## psychonetrics model comparison ##\n",
  #          "\t\t####################################\n\n"))
  # psychonetrics_print_logo()
  # No awesome header :(
  
  print.data.frame(x,row.names=FALSE)
  
  cat("\nNote: Chi-square difference test assumes models are nested.")
}
