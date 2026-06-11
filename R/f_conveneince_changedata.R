changedata <- function(x, data, covs, nobs, means, groups, missing = "listwise"){
  if (x@distribution != "Gaussian"){
    stop("Only Gaussian distirbution models supported.")
  }
  if (x@model == "var1"){
    stop("VAR models not supported yet.")
  }
  
  
 # Set computed to false:
  x@computed <- FALSE
  
  # Check if labels are the same (only possible when raw data is supplied;
  # the covs/means/nobs path is validated downstream by samplestats):
  if (!missing(data)){
    if (!all(x@sample@variables$label %in% colnames(data))){
      stop("New dataset does not contain the same variables as used in model.")
    }
  }

  # Variables:
  vars <-  x@sample@variables$label

  # Estimator:
  estimator <- x@estimator

  # Collect the data-source arguments: either raw `data`, or summary
  # statistics (covs/means/nobs), mirroring how the model constructors accept
  # either form of input. Only the supplied arguments are forwarded so the
  # constructors' own missing() handling kicks in for the rest.
  dataArgs <- list()
  if (!missing(data))   dataArgs$data   <- data
  if (!missing(groups)) dataArgs$groups <- groups
  if (!missing(covs))   dataArgs$covs   <- covs
  if (!missing(means))  dataArgs$means  <- means
  if (!missing(nobs))   dataArgs$nobs   <- nobs

  # Data set:
  x@sample <- do.call(samplestats, c(dataArgs, list(
              vars = vars,
              missing  = ifelse(estimator == "FIML","pairwise",missing),
              fimldata = estimator == "FIML")))

  # Number of observations (degrees-of-freedom base). Mirror the model
  # constructor: variances/covariances per group, plus the mean block only when
  # the model has a mean structure (otherwise df would be corrupted for
  # meanstructure = FALSE models):
  nVar <- nrow(x@sample@variables)
  nGroup <- nrow(x@sample@groups)
  x@sample@nobs <-
    nVar * (nVar+1) / 2 * nGroup +              # Variances/covariances per group
    x@meanstructure * nVar * nGroup             # Means per group (if present)
  
  # Add baseline and saturated:
  
  # Form baseline model:
  # NOTE: equal = x@equal forwards the legacy (deprecated) @equal slot (a
  # character vector of matrix names). The baseline is always a varcov() model,
  # so only equality constraints on matrices that exist in varcov (sigma, omega,
  # kappa, lowertri, ...) carry over; constraints on framework-specific matrices
  # (e.g. lambda in lvm) are silently dropped. Deriving equal= from the canonical
  # parameters table would be more robust but non-trivial; documented limitation.
  x@baseline_saturated$baseline <- do.call(varcov, c(dataArgs, list(
                                              type = "chol",
                                              lowertri = "diag",
                                              vars = vars,
                                              missing = missing,
                                              equal = x@equal,
                                              estimator = estimator,
                                              baseline_saturated = FALSE)))
  
  # Add model:
  # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
  
  
  ### Saturated model ###
  # No `equal = x@equal`: the saturated reference is unconstrained per
  # group; cross-group equality belongs to the target/baseline only.
  x@baseline_saturated$saturated <- do.call(varcov, c(dataArgs, list(
                                               type = "chol",
                                               lowertri = "full",
                                               vars = vars,
                                               missing = missing,
                                               estimator = estimator,
                                               baseline_saturated = FALSE)))
  
  # if not FIML, Treat as computed:
  if (estimator != "FIML"){
    x@baseline_saturated$saturated@computed <- TRUE
    
    # FIXME: TODO
    x@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$saturated),x@baseline_saturated$saturated)      
  }
  
  return(x)
}
