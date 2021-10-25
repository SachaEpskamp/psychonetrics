labtoind <- function(x,row,col,Mat,symmetrical=FALSE){
  symmetrical <- x@matrices$symmetrical[x@matrices$name == Mat]

  # Obtain all labels in rows:
  if (symmetrical){
    vars <- rbind(
      x@parameters %>% filter(.data[['matrix']] == Mat) %>% select(ind = .data[['row']], var = .data[['var1']]),
      x@parameters %>% filter(.data[['matrix']] == Mat) %>% select(ind = .data[['col']], var = .data[['var2']])
    )
    vars <- vars[!duplicated(vars),]
    
    row <- vars$ind[match(row,vars$var)]
    col <- vars$ind[match(col,vars$var)]
  } else{
    
    rows <- x@parameters %>% filter(.data[['matrix']] == Mat) %>% group_by(.data[["row"]]) %>%
      summarize(ind = unique(.data[['row']]), var = unique(.data[['var1']]))
    cols <-  x@parameters %>% filter(.data[['matrix']] == Mat) %>% group_by(.data[["col"]]) %>%
      summarize(ind = unique(.data[['col']]), var = unique(.data[['var2']]))
    
    
    row <- rows$ind[match(row,rows$var)]
    col <- cols$ind[match(col,cols$var)]
  }
 
  return(list(
    row=row,col=col
  ))
}