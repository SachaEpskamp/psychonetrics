labtoind <- function(x,row,col,Mat,symmetrical=FALSE){
  symmetrical <- x@matrices$symmetrical[x@matrices$name == Mat]

  # Obtain all labels in rows:
  if (symmetrical){
    vars <- rbind(
      x@parameters %>% filter(.data[['matrix']] == Mat) %>% select(ind = .data[['row']], var = .data[['var1']]),
      x@parameters %>% filter(.data[['matrix']] == Mat) %>% select(ind = .data[['col']], var = .data[['var2']])
    )
    vars <- vars[!duplicated(vars),]
    
    if (is.character(row)){
      unmatched <- row[!is.na(row) & !(row %in% vars$var)]
      if (length(unmatched) > 0){
        stop(paste0("The following variable label(s) do not occur in matrix '",Mat,"': ", paste0(unique(unmatched), collapse=", ")))
      }
    }
    if (is.character(col)){
      unmatched <- col[!is.na(col) & !(col %in% vars$var)]
      if (length(unmatched) > 0){
        stop(paste0("The following variable label(s) do not occur in matrix '",Mat,"': ", paste0(unique(unmatched), collapse=", ")))
      }
    }

    row <- vars$ind[match(row,vars$var)]
    col <- vars$ind[match(col,vars$var)]
  } else{
    
    rows <- x@parameters %>% filter(.data[['matrix']] == Mat) %>% group_by(.data[["row"]]) %>%
      summarize(ind = unique(.data[['row']]), var = unique(.data[['var1']]))
    cols <-  x@parameters %>% filter(.data[['matrix']] == Mat) %>% group_by(.data[["col"]]) %>%
      summarize(ind = unique(.data[['col']]), var = unique(.data[['var2']]))

    if (is.character(row)){
      unmatched <- row[!is.na(row) & !(row %in% rows$var)]
      if (length(unmatched) > 0){
        stop(paste0("The following variable label(s) do not occur in the rows of matrix '",Mat,"': ", paste0(unique(unmatched), collapse=", ")))
      }
    }
    if (is.character(col)){
      unmatched <- col[!is.na(col) & !(col %in% cols$var)]
      if (length(unmatched) > 0){
        stop(paste0("The following variable label(s) do not occur in the columns of matrix '",Mat,"': ", paste0(unique(unmatched), collapse=", ")))
      }
    }

    row <- rows$ind[match(row,rows$var)]
    col <- cols$ind[match(col,cols$var)]
  }
 
  return(list(
    row=row,col=col
  ))
}