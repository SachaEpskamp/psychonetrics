labtoind <- function(x,row,col,Mat,symmetrical=FALSE){
  symmetrical <- x@matrices$symmetrical[x@matrices$name == Mat]

  # Obtain all labels in rows:
  if (symmetrical){
    vars <- rbind(
      x@parameters %>% filter_(~matrix == Mat) %>% select_(ind = ~row, var = ~var1),
      x@parameters %>% filter_(~matrix == Mat) %>% select_(ind = ~col, var = ~var2)
    )
    vars <- vars[!duplicated(vars),]
    
    row <- vars$ind[match(row,vars$var)]
    col <- vars$ind[match(col,vars$var)]
  } else{
    
    rows <- x@parameters %>% filter_(~matrix == Mat) %>% group_by_("row") %>%
      summarize_(ind = ~unique(row), var = ~unique(var1))
    cols <-  x@parameters %>% filter_(~matrix == Mat) %>% group_by_("col") %>%
      summarize_(ind = ~unique(col), var = ~unique(var2))
    
    
    row <- rows$ind[match(row,rows$var)]
    col <- cols$ind[match(col,cols$var)]
  }
 
  return(list(
    row=row,col=col
  ))
}