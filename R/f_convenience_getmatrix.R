# obtain a model matrix:
getmatrix <- function(x,matrix,group){
  if (!is(x,"psychonetrics")){
    stop("Input must be a 'psychonetrics' object.")
  }
  
  # check matrix arg:
  if (missing(matrix)){
    stop("'matrix' argument may not be missing.")
  }
  if (!is.character(matrix) || length(matrix) > 1 || !matrix %in% names(x@modelmatrices[[1]])){
    stop("'matrix' argument is not a character string of a single matrix in the model.")
  }
  
  # If group is missing, all groups:
  if (missing(group)){
    group <- x@sample@groups$label
  }
  
  # If group is number, get name:
  if (is.numeric(group)){
    group <- x@sample@groups$label[match(group,x@sample@groups$id)]
  }
  
  # Form group ID:
  groupID <- x@sample@groups$id[match(group,x@sample@groups$label)]
  
  # Obtain matrices:
  mats <- lapply(x@modelmatrices[groupID],function(x)as.matrix(x[[matrix]]))
  names(mats) <- group
  
  # If length = 1, only return the matrix:
  if (length(mats)==1){
    mats <- mats[[1]]
  }
  return(mats)
}