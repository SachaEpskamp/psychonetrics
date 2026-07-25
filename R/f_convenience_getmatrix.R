# obtain a model matrix:
getmatrix <- function(x,matrix,group,threshold=FALSE,
                      alpha = 0.01, 
                      adjust = c( "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                      mode = c("tested","all"),
                      diag = TRUE,
                      time
                      # inclusion = c("any","positive","negative")
                      ){
  mode <- match.arg(mode)


  
  # # Check input:
  # if (is(x,"psychonetrics_bootstrap")){
  #   boot_aggregate <- TRUE
  #   if (!(is.logical(threshold) && threshold == FALSE)){
  #     stop("Thresholding not supported for bootstrap aggregates")
  #   }
  # } else if (!is(x,"psychonetrics")){
  #   stop("Input must be a 'psychonetrics' object.")
  # } else {
  #   boot_aggregate <- FALSE
  # }
  
  if (!is(x,"psychonetrics")){
    stop("Input must be a 'psychonetrics' object.")
  }

  # --- RI-CLPM / RI-CLPN block requests ---
  # ri_clpm() is built on the lvm framework, so its matrices span all waves at
  # once. Requests such as getmatrix(x, "omega_zeta_within", time = 2) or
  # getmatrix(x, "omega_zeta_between") are translated here into (real matrix,
  # rows, cols); the real matrix then runs through the machinery below exactly
  # as before and the block is cut out at the very end, so that thresholding
  # and multiple-comparison adjustment behave identically to extracting the
  # block from the full matrix by hand. Plain names without 'time' keep their
  # old meaning (the full matrix), so this is backward compatible.
  riclpm_spec <- NULL
  if (identical(x@submodel, "RI_CLPM") && !missing(matrix) &&
      is.character(matrix) && length(matrix) == 1){
    riclpm_spec <- if (missing(time)) riclpm_resolve(x, matrix) else riclpm_resolve(x, matrix, time)
    if (!is.null(riclpm_spec)) matrix <- riclpm_spec$real
  } else if (!missing(time)){
    stop("The 'time' argument is only available for ri_clpm() / ri_clpn() models.")
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

  # --- RI-CLPM: partial directed correlations of a transition ---
  # The RI-CLPM has no PDC model matrix (PDC is not one of its
  # parameterizations), so it is derived per transition from the temporal
  # block and the innovation precision of the target wave. Thresholding
  # follows the same rule as elsewhere: the PDC of a temporal effect that is
  # not significant is set to zero.
  if (!is.null(riclpm_spec) && riclpm_spec$kind == "PDC"){
    out <- lapply(riclpm_spec$views, function(v){
      pdc <- riclpm_PDC(x, v, groupID)
      names(pdc) <- group
      if (!identical(threshold, FALSE)){
        b <- getmatrix(x, "beta", group, threshold = threshold, alpha = alpha,
                       adjust = adjust, mode = mode, diag = TRUE, time = v$waves[1])
        if (!is.list(b)) b <- list(b)
        for (g in seq_along(pdc)) pdc[[g]] <- pdc[[g]] * t(b[[g]] != 0)
      }
      if (!diag) for (g in seq_along(pdc)) diag(pdc[[g]]) <- 0
      if (length(pdc) == 1) pdc[[1]] else pdc
    })
    names(out) <- vapply(riclpm_spec$views, "[[", character(1), "label")
    if (length(out) == 1) return(out[[1]])
    return(out)
  }

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
  
  # --- RI-CLPM: cut the requested block(s) out of the full matrix ---
  # Done after thresholding so that a block is exactly the corresponding
  # submatrix of the thresholded full matrix (identical p-value adjustment),
  # but before 'diag' so that diag = FALSE removes the diagonal OF THE BLOCK
  # (for a temporal block that is the autoregression, which does not lie on
  # the diagonal of the full beta matrix).
  if (!is.null(riclpm_spec)){
    blocks <- lapply(riclpm_spec$views, function(v){
      bl <- lapply(mats, riclpm_subset, view = v, spec = riclpm_spec)
      if (!diag){
        bl <- lapply(bl, function(m){
          if (nrow(m) == ncol(m)) diag(m) <- 0
          m
        })
      }
      if (length(bl) == 1) bl[[1]] else bl
    })
    names(blocks) <- vapply(riclpm_spec$views, "[[", character(1), "label")
    if (length(blocks) == 1) return(blocks[[1]])
    return(blocks)
  }

  # Remove diagonal if needed:
  if (!diag){
    for (g in seq_along(mats)){
      if (ncol(mats[[g]]) == nrow(mats[[g]]))
      diag(mats[[g]]) <- 0
    }
  }

  # If length = 1, only return the matrix:
  if (length(mats)==1){
    mats <- mats[[1]]
  }


  return(mats)
}