# Generate multiple par tables:
generateAllParameterTables <- function(...){
  curMaxPar <- 0
  dots <- list(...)
  res <- list()
  for (i in seq_along(dots)){
    res[[i]] <- do.call(generateParameterTable,c(dots[[i]],list(curMaxPar=curMaxPar)))
    curMaxPar <- max(res[[i]]$par)
  }
  do.call(rbind,res)
}

generateParameterTable <- function(x, mat, op, curMaxPar, symmetrical = FALSE, sampletable, rownames, colnames, rowid, colid){
  # rowid and colid can be missing:
  if (missing(rowid)){
    rowid <- seq_along(rownames)
  }
  if (missing(colid)){
    colid <- seq_along(colnames)
  }
  
  # Indices to use:
  if (symmetrical){
    ind <- lower.tri(x,diag=TRUE)
  } else {
    ind <- rep(TRUE,length(x))
  }
  # Obtain row:
  row <- slice.index(x,1)[ind]
  
  if (length(dim(x)) == 3){
    # Obtain col:
    col <- slice.index(x,2)[ind]
    # Obtain group:
    group <- slice.index(x,3)[ind]    
  } else {
    # Obtain col:
    col <- NA
    # Obtain group:
    group <- slice.index(x,2)[ind]
  }

  # Est = 0 if matrix != symmetrical, 1 for diagonals if matrix = symmetrical:
  if (symmetrical){
    est <- array(diag(nrow(x)),dim(x))[ind]
  } else {
    est <- 0
  }
  
  # Fixed are values that are zero:
  fixed <- x[ind] == 0
  
  # Fixed values are always set to zero:
  est[fixed] <- 0
  
  # Initiate pars at zero:
  par <- rep(0,length=length(est))
  
  # All ones, make increasing sequence:
  par[x[ind]==1] <- curMaxPar + seq_len(sum(x[ind]==1))
  curMaxPar <- max(par)
  
  # Now loop if needed:
  if (any(x[ind] > 1)){
    for (p in unique(x[ind & x[ind]!=0 & x[ind]!=1])){
      par[x[ind] == p] <- curMaxPar + 1
      curMaxPar <- curMaxPar + 1
    }
  }
   
  # var1:
  var1 <- rownames[row]
  var1_id <- rowid[row]
  if (length(dim(x))>2){
    var2 <- rownames[col]
    var2_id <- colid[row]    
  } else {
    var2 <- NA
    var2_id <- NA
  }
  group_name <- sampletable@groups$label[group]
  group_id <- group

  # Make the table:
  table <- data.frame(
    var1 = var1,
    var1_id = var1_id,
    op = op,
    var2 = var2,
    var2_id = var2_id,
    est = est,
    std = NA,
    se = NA,
    p = NA,
    matrix = mat,
    row = row,
    col = col,
    par = par,
    group = group_name,
    group_id = group_id,
    fixed = fixed
  )
  
  return(table)
}