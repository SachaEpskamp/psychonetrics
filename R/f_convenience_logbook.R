# psychonetrics print method
 
setMethod(f = "show",
 signature = "psychonetrics_log",
definition = function(object){
  
  time <- object@time
  event <- object@event
  
  cat (paste0(
    "Event logged at ", format(time, "%Y-%m-%d %H:%M"),":\n\t",event,"\n")
  )
})
  

# Logbook:
print.psychonetrics_log <- function(x,...){
  
  # Extract the POSIXct timestamps (using do.call(c,...) preserves the
  # POSIXct class, unlike sapply which would strip it to a bare numeric):
  times <- do.call(c, lapply(x, slot, "time"))

  df <- data.frame(
    time =  format(times, "%Y-%m-%d %H:%M"),
   events =  sapply(x, 'slot', "event")
  )
  
  # Print:
  for (i in seq_len(nrow(df))){
    cat(paste0("Event logged at ",df$time[i],": \n\t",df$events[i],"\n\n"))
  }
  
  return(df)
}

# Convenience function:
logbook <- function(x,log=TRUE){
  stopifnot(is(x,"psychonetrics"))
  
  # Add final entry for when the logbook is accessed:
  if (log){
    x <- addLog(x,"Logbook retrieved")  
  }
  
  # Return logook:
  return(x@log)
}