# Generate multiple par tables:
generateAllParameterTables <- function(...){
  curMaxPar <- 0
  dots <- list(...)
  res <- list()
  for (i in seq_along(dots)){
    res[[i]] <- do.call(generateParameterTable,c(dots[[i]],list(curMaxPar=curMaxPar)))
    curMaxPar <- max(res[[i]]$partable$par)
  }
  res <- list(
    partable = do.call(rbind,lapply(res,"[[","partable")),
    mattable = do.call(rbind,lapply(res,"[[","mattable"))
  )
  
  # Order parameter table by group id:
  res$partable <- res$partable %>% arrange_(~group_id)
  # Relabel the parameter labels to be a bit more consistent:
  res$partable$par[res$partable$par!=0] <- as.numeric(factor(res$partable$par[res$partable$par!=0], 
            levels = unique(res$partable$par[res$partable$par!=0])))
  
  # Return:
  res
}

generateParameterTable <- function(x, mat, op, curMaxPar, symmetrical = FALSE, sampletable, rownames, colnames, rowid, colid, sparse = FALSE, posdef = FALSE, diag0=FALSE, diagonal = FALSE){
  # rowid and colid can be missing:
  if (missing(rowid)){
    rowid <- seq_along(rownames)
  }
  if (missing(colid)){
    colid <- seq_along(colnames)
  }
  
  # Indices to use:
  if (diagonal){
    ind <- array(FALSE,dim=dim(x))
    for (i in 1:dim(ind)[3]){
      ind[,,i] <- diag(nrow(ind[,,i])) == 1
    }
  } else if (symmetrical){
    ind <- array(TRUE,dim=dim(x))
    for (i in 1:dim(ind)[3]){
      ind[,,i] <- lower.tri(ind[,,i],diag=!diag0)
    }
    # ind <- lower.tri(x,diag=TRUE)
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
    col <- 1
    # Obtain group:
    group <- slice.index(x,2)[ind]
  }

  # Est = 0 if matrix != symmetrical, 1 for diagonals if matrix = symmetrical:
  if (symmetrical){
    est <- array(diag(nrow(x)),dim(x))[ind]
  } else {
    est <- rep(0,length(x))
  }
  
  # Fixed are values that are zero:
  fixed <- x[ind] == 0
  
  # Fixed values are always set to zero:
  est[fixed] <- 0
  
  # Initiate pars at zero:
  par <- rep(0,length=length(est))
  
  # All ones, make increasing sequence:
  if (any(x[ind]==1)){
    par[x[ind]==1] <- curMaxPar + seq_len(sum(x[ind]==1))
    curMaxPar <- max(par)    
  }


  # Now loop if needed:
  if (any(x[ind] > 1)){
    for (p in unique(x[ind & x!=0 & x!=1])){
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
  partable <- data.frame(
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
    fixed = fixed,
    mi = NA, # Modification index
    pmi = NA, #p-value modification index
    mi_equal = NA, # Modification index constraning groups to be equal
    pmi_equal = NA, #p-value modification index constraining groups to be equal
    minimum = -Inf,
    maximum = Inf,
    stringsAsFactors = FALSE
  )
  
  #FIXME: Temporary fix to make diagonals of symmetrical matrices non-negative:
  if (symmetrical){
    partable$minimum[partable$var1_id == partable$var2_id] <- 1e-14
  }
  
  # Table for matrices:
  mattable <- data.frame(
    name = mat,
    nrow = dim(x)[1],
    ncol = ifelse(length(dim(x)) > 2,dim(x)[2],1),
    ngroup = max(group_id),
    symmetrical = symmetrical,
    sparse = sparse,
    posdef=posdef,
    diagonal= diagonal
  )
  
  return(list(
    partable = partable,
    mattable = mattable
  ))
}