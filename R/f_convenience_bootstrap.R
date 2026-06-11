# Bootstrap the data:
bootstrap <- function(x, replacement = TRUE, proportion = 1, verbose = TRUE, storedata = FALSE, baseline_saturated = TRUE){
  if (x@distribution != "Gaussian"){
    stop("Only Gaussian distirbution models supported.")
  }
  
  # Set computed to false:
  x@computed <- FALSE
  
  # Obtain the raw data:
  if (is.null(x@sample@rawdata) || nrow(x@sample@rawdata) == 0){
    stop("Raw data is needed. Please run input with 'storedata = TRUE'")
  }
  data <- x@sample@rawdata

  # Data attributes (variable names, grouping column, missingness handling).
  # Capture them up front: row resampling and rbind() below may not carry
  # custom attributes reliably across R versions.
  varsAttr    <- attr(data, "vars")
  groupVar    <- attr(data, "groups")
  missingAttr <- attr(data, "missing")

  # bootstrap the data:
  # Resample STRATIFIED within each group, preserving each group's size, and
  # rebuild the data in the ORIGINAL group order (as stored in
  # x@sample@groups$label). Pooled resampling would (a) let group sizes drift
  # between replications and (b) make samplestats() derive the group order from
  # whichever group the first resampled row happens to belong to -- which can
  # permute the group ids relative to x@parameters and fit group 1's parameters
  # to another group's data.
  groupLabels <- x@sample@groups$label
  origData <- data
  data <- do.call(rbind, lapply(groupLabels, function(g){
    rows <- which(origData[[groupVar]] == g)
    take <- round(proportion * length(rows))
    origData[sample(rows, take, replace = replacement), , drop = FALSE]
  }))

  # New data:
  x@sample <- samplestats(
    data = data, # Dataset
    vars = varsAttr, # character indicating the variables Extracted if missing from data - group variable
    groups = groupVar, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
    # covs, # alternative covs (array nvar * nvar * ngroup)
    # means, # alternative means (matrix nvar * ngroup)
    # nobs, # Alternative if data is missing (length ngroup)
    missing =  missingAttr,
    fimldata = x@estimator == "FIML",
    verbose = verbose,
    storedata = storedata
  )

  # Number of observations (degrees-of-freedom base). Mirror the model
  # constructor: variances/covariances per group, plus the mean block only when
  # the model has a mean structure (otherwise df would be corrupted for
  # meanstructure = FALSE models):
  nVar <- nrow(x@sample@variables)
  nGroup <- nrow(x@sample@groups)
  x@sample@nobs <-
    nVar * (nVar+1) / 2 * nGroup +              # Variances/covariances per group
    x@meanstructure * nVar * nGroup             # Means per group (if present)
  
  # Reset start values:
  x <- emergencystart(x)
  
  # Recompute baseline and saturated:
  
  # Form baseline model:
  if (baseline_saturated){
    x@baseline_saturated$baseline <- varcov(data,
                                            type = "chol",
                                            lowertri = "diag",
                                            vars = varsAttr,
                                            groups = groupVar,
                                            missing = missingAttr,
                                            equal = x@equal,
                                            estimator = x@estimator,
                                            baseline_saturated = FALSE)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    # No `equal = x@equal`: the saturated reference is unconstrained per
    # group; cross-group equality belongs to the target/baseline only.
    x@baseline_saturated$saturated <- varcov(data,
                                             type = "chol",
                                             lowertri = "full",
                                             vars = varsAttr,
                                             groups = groupVar,
                                             missing = missingAttr,
                                             estimator = x@estimator,
                                             baseline_saturated = FALSE)
    
    # if not FIML, Treat as computed:
    if (x@estimator != "FIML"){
      x@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      x@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$saturated),x@baseline_saturated$saturated)      
    }
  } else {
    x@baseline_saturated <- list()
  }
  
  
  # Return:
  return(x)
}