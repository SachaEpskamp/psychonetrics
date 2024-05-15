# obtain a model matrix:
getmatrix <- function(x,matrix,group,threshold=FALSE,
                      alpha = 0.01, 
                      adjust = c( "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                      mode = c("tested","all"),
                      diag = TRUE){
  mode <- match.arg(mode)
  

  # Check input:
  if (!is(x,"psychonetrics")){
    stop("Input must be a 'psychonetrics' object.")
  }
  
  # IF matrix is PDC, run recursive:
  if (matrix == "PDC" & !identical(threshold,FALSE)){
    beta <- getmatrix(x,"beta",group,threshold=threshold,alpha=alpha,adjust=adjust,mode=mode,diag=diag)
    PDC <- getmatrix(x,"PDC",group,threshold=FALSE,mode=mode,diag=diag)
    
    if (is.list(PDC)){
      reslist <- PDC
      for (g in seq_along(reslist)){
        reslist[[g]] <- PDC[[g]] * t(beta[[g]]!=0)
      }
    } else {
      reslist <- PDC * t(beta!=0)
    }
    return(reslist)
  }
  
  # Extract verbose:
  verbose <- x@verbose
  
  
  # If not run, run model:
  if (!x@computed){
    warning("Model has not been computed! Returning start-values.")
  #   x <- x %>% runmodel(verbose = verbose)
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
  } else   # If group is number, get name:
    if (is.numeric(group)){
      group <- x@sample@groups$label[match(group,x@sample@groups$id)]
    }
  
  
  
  # Form group ID:
  groupID <- x@sample@groups$id[match(group,x@sample@groups$label)]
  
  # Obtain matrices:
  mats <- lapply(x@modelmatrices[groupID],function(x)as.matrix(x[[matrix]]))
  names(mats) <- group
  
  
  # threshold:
  if (isTRUE(threshold)){
    reps <- 1000
    nCores <- 1
    bootstrap <- FALSE
    adjust <- match.arg(adjust)
    
    # Check if the matrix is modeled:
    # matrix <- matrix
    # if (matrix == "PDC"){
    #   # FIXME: 
    #   stop("PDC thresholding is not directly implemented. For now, use getmatrix(mod,'PDC') * t(getmatrix(mod,'beta',threshold=TRUE,...)!=0)")
    #   # matrix <- "beta"
    # }
    if (!matrix %in% x@matrices$name){
      stop("Matrix is not modeled and can therefore not be thresholded.")
    }

    # obtain p-values:
    pValues <- adjust_p_values(x,
                               alpha = alpha, 
                               adjust = adjust,
                               matrices = matrix,
                               mode = mode,
                               reps = reps,
                               nCores = nCores,
                               bootstrap = bootstrap,
                               verbose = verbose)

    sig <- !is.na(pValues) & pValues <= alpha # & (seq_len(nrow(x@parameters)) %in% whichTest)
    
    # Threshold all nonsig to zero:
    copy_pars <- x@parameters
    copy_pars$est_thresholded <- sig * copy_pars$est
    
    # Loop over groups:
    for (g in seq_along(mats)){
      # Obtain the relevant parameter table:
      subPars <- copy_pars[copy_pars$matrix == matrix & copy_pars$group_id == groupID[g], ]
      
      # Loop over the parameters and overwrite:
      for (i in seq_len(nrow(subPars))){
        mats[[g]][subPars$row[i],subPars$col[i]] <- subPars$est_thresholded[i]
        if (x@matrices$symmetrical[x@matrices$name == matrix]){
          mats[[g]][subPars$col[i],subPars$row[i]] <- subPars$est_thresholded[i]
        }
      }
    }
  } else if (is.numeric(threshold)){
    for (g in seq_along(mats)){
      mats[[g]][abs(mats[[g]]) < threshold] <- 0
    }
  }
  
  # Remove diagonal if needed:
  if (!diag){
    for (g in seq_along(mats)){
      if (ncol(mats[[g]]) == nrow(mats[[g]]))
      mats[[g]][diag(mats[[g]])] <- 0
    }
  }
  
  # If length = 1, only return the matrix:
  if (length(mats)==1){
    mats <- mats[[1]]
  }
  
  
  return(mats)
}