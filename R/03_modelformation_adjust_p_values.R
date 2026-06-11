adjust_p_values <- function(x,
          alpha = 0.01, 
          adjust = c( "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
          matrices,
          mode = c("tested","all"),
          reps = 1000,
          nCores = 1,
          bootstrap = FALSE,
          verbose){
  

  # input:
  stopifnot(!missing(matrices))
  mode <- match.arg(mode)
  adjust <- match.arg(adjust)
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  # Bootstrap (currently not used):
  if (bootstrap && all(is.na(x@parameters$boot_p))){
    if (verbose) message("Bootstrapping SEs...")
    x <- x %>% bootstrap_SEs(nCores = nCores, reps = reps,verbose = verbose)
  }
  
  # Whcih cols?
  if (bootstrap){
    secol <- "se_boot"
    pcol <- "p_boot" 
  } else {
    secol <- "se"
    pcol <- "p"
  }
  
  
  # Which parameters to test:
  # FIXME: Not sure why original version removes diagonal elements?
  # whichTest <- which(x@parameters$matrix %in% matrices & !x@parameters$fixed & (NAtoTRUE(x@parameters$var1_id!=x@parameters$var2_id) | x@parameters$matrix == "beta"))
  whichTest <- which(x@parameters$matrix %in% matrices & !x@parameters$fixed & (NAtoTRUE(x@parameters$var1_id!=x@parameters$var2_id) | !grepl("omega",x@parameters$matrix)))
  
  # Number of tests:
  nTest <- length(unique(x@parameters$par[whichTest]))
  
  # If no tests, return NA vector:
  if (nTest == 0){
    return(rep(NA, nrow(x@parameters)))
  }
  
  # Parameter indices (equality-constrained parameters share a par number > 0
  # across rows/groups and must be counted only once in the multiple-comparison
  # correction; fixed parameters have par == 0 and are not tested here):
  par <- x@parameters$par

  if (mode == "all"){
    pValues <- rep(NA, nrow(x@parameters))
    # Adjust over unique (free) parameters, then map back to all rows:
    rows <- which(par > 0)
    if (length(rows) > 0){
      u <- !duplicated(par[rows])
      padj_u <- p.adjust(x@parameters[[pcol]][rows][u], method = adjust)
      # Map each row to its parameter's adjusted p-value:
      pValues[rows] <- padj_u[match(par[rows], par[rows][u])]
    }
  } else {
    pValues <- rep(NA, nrow(x@parameters))
    # Adjust over unique tested parameters, then map back to all tested rows:
    u <- !duplicated(par[whichTest])
    padj_u <- p.adjust(x@parameters[[pcol]][whichTest][u], method = adjust)
    pValues[whichTest] <- padj_u[match(par[whichTest], par[whichTest][u])]
  }

  # Return:
  pValues
}