# This function takes a bootnet object, and extracts the data and network:
frombootnet <- function(x, model, ...){
  # Extract network:
  if (is(x,"bootnet")){
    x <- x$sample
  }
  
  # If not the right class, stop:
  if (!is(x,"bootnetResult")){
    stop("Input object must be of class 'bootnetResult'")
  }

  # If model is missing, set:
  if (missing(model)){
    if (x$default %in% c("EBICglasso", "pcor", "huge", "adalasso", "ggmModSelect")){
      model <- "ggm"
    } else if (x$default %in% c("IsingSampler","IsingFit")){
      model <- "Ising"
    } else {
      stop("'model' could not be determined")
    }
  }
  
  # FIXME: Currently only ggm supported
  if (model != "ggm"){
    stop("Only model = 'ggm' is supported")
  }
  
  # Make model:
  if (model == "ggm"){
    mod <- ggm(x$data, omega = 1*(x$graph!=0),vars = x$labels, ...)
  }
  
  # Add log:
  mod <- addLog(mod, paste0("Model created from bootnet object! Default used: ",x$default)) 
  # Remove first entry (created):
  mod@log[[1]] <- NULL
  
  return(mod)
}