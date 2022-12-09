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
  
  # If no tests, break:
  if (nTest == 0){
    return(x)
  }
  
  # old method:
  if (mode == "all"){
    pValues <- p.adjust(x@parameters[[pcol]],method = adjust) 
  } else {
    pValues <- rep(NA,nrow(x@parameters))
    pValues[whichTest] <- p.adjust(x@parameters[[pcol]][whichTest],method = adjust) 
  }
  
  # Return:
  pValues
}