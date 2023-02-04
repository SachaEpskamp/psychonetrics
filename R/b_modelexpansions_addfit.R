# Add fit measures to psychonetrics object!



# Computes fit measures
addfit <- function(
 x, #, ebicTuning = 0.25
 verbose
){
  if (missing(verbose)){
    verbose <- x@verbose
  }
  
  if (verbose){
    message("Adding fit measures...")
  }
  
  # If not computed, stop:
  if (!x@computed){
    stop("Model has not yet been run. Use runmodel(object) first!")
  }


  # Sample size:
  sampleSize <- sum(x@sample@groups$nobs)
  
  # Fitmeasures list:
  fitMeasures <- list()
  
  
  # log likelihoods:
  # Saturated:
  if (x@estimator %in% c("FIML","ML")){
    if (!is.null(x@baseline_saturated$saturated)){
      satLL <- psychonetrics_logLikelihood(x@baseline_saturated$saturated)    
    } else {
      satLL <- NA
    }
    # Baseline:
    if (!is.null(x@baseline_saturated$baseline)){
      basLL <- psychonetrics_logLikelihood(x@baseline_saturated$baseline)
    } else {
      basLL <- NA
    }
    
    # Model:
    LL <-  psychonetrics_logLikelihood(x)
  } else {
    satLL <- NA
    basLL <- NA
    LL <- NA
  }

  # Add to list:
  fitMeasures$logl <- LL
  fitMeasures$unrestricted.logl <- satLL
  fitMeasures$baseline.logl <- basLL
  
  # Number of variables:
  fitMeasures$nvar <- nVar <- nrow(x@sample@variables)
  
  # Number of observations:
  fitMeasures$nobs <- x@sample@nobs
    
  # Number of parameters:
  fitMeasures$npar <- max(x@parameters$par)
  
  # Degrees of freedom:
  fitMeasures$df <- fitMeasures$nobs - fitMeasures$npar    
  if (!is.null(x@baseline_saturated$saturated)){
    fitMeasures$df <- fitMeasures$df  - (x@baseline_saturated$saturated@sample@nobs - max(x@baseline_saturated$saturated@parameters$par))
  } 

  # Compute Fmin:
  fitMeasures$objective <- x@objective
  
  # Likelihood ratio:
  if (x@estimator %in% c("FIML","ML")){
    fitMeasures$chisq <- -2 * (LL - satLL)
  } else  if (x@estimator %in% c("WLS","DWLS","ULS")){
    fitMeasures$chisq <- x@objective  * (sampleSize)
  } 
  fitMeasures$pvalue <- pchisq(fitMeasures$chisq, fitMeasures$df, lower.tail = FALSE)
  
  # Some pars:
  Tm <- fitMeasures$chisq
  dfm <- fitMeasures$df
  
  # Baseline model:
  if (!is.null(x@baseline_saturated$baseline) && x@baseline_saturated$baseline@computed){
    if (length(x@baseline_saturated$baseline@objective) == 0){
      x@baseline_saturated$baseline@objective <- psychonetrics_fitfunction(parVector(x@baseline_saturated$baseline),x@baseline_saturated$baseline)
    }
    
    # fitMeasures$fmin_baseline <- x@baseline_saturated$baseline@objective
    # fitMeasures$baseline.chisq <-  sampleSize * fitMeasures$fmin_baseline
    if (x@estimator%in% c("FIML","ML")){
      fitMeasures$baseline.chisq <-  -2 * (basLL - satLL)
    } else  if (x@estimator %in% c("WLS","DWLS","ULS")){
      fitMeasures$baseline.chisq <- x@baseline_saturated$baseline@objective  * (sampleSize)
    } 

    fitMeasures$baseline.df <- fitMeasures$nobs - max(x@baseline_saturated$baseline@parameters$par)
    fitMeasures$baseline.pvalue <- pchisq(fitMeasures$baseline.chisq, fitMeasures$baseline.df, lower.tail = FALSE)
    
    # Incremental Fit Indices
    Tb <- fitMeasures$baseline.chisq
    
    dfb <- fitMeasures$baseline.df
    # 
    # t1 <- (X2 - df)*df.null
    # t2 <- (X2.null - df.null)*df
    # if(df > 0 && abs(t2) > 0) {
    #   indices["tli"] <- indices["nnfi"] <- 1 - t1/t2
    # } else {
    #   indices["tli"] <- indices["nnfi"] <- 1
    # }
    
    fitMeasures$nfi <- (Tb - Tm) / Tb

    # Stop here if baseline is not good:
    if (is.null(dfb) || !is.finite(dfb) || !is.finite(Tb)){
      return(x)
    }
    
    if(dfb > 0 && Tb > 0) {
      t1 <- Tb - Tm
      t2 <- Tb
      fitMeasures$pnfi <- (dfm/dfb) * t1/t2
    } else {
      fitMeasures$pnfi <- as.numeric(NA)
    }
    
    fitMeasures$tli <-  (Tb/dfb - Tm/dfm) / (Tb/dfb - 1) 

    
    t1 <- (Tm - dfm)*dfb
    t2 <- (Tb - dfb)*dfm
    fitMeasures$nnfi <- ifelse(dfm > 0 & abs(t2) > 0, 1 - t1/t2, 1)
    

    
    fitMeasures$rfi <-  (Tb/dfb - Tm/dfm) / (Tb/dfb ) 
    fitMeasures$ifi <-  (Tb - Tm) / (Tb - dfm)
    fitMeasures$rni <-  ((Tb- dfb) - (Tm - dfm)) / (Tb - dfb)
    fitMeasures$cfi <- ifelse(dfm > Tm, 1, 1 - (Tm - dfm)/(Tb - dfb))
    
    } else {
    warning("No baseline model found, cannot add incremental fit indices...")
      fitMeasures$fmin_baseline <- NA
      fitMeasures$baseline.chisq <- NA
      fitMeasures$baseline.df <- NA
      fitMeasures$baseline.pvalue <- NA
      
      # Incremental Fit Indices
      fitMeasures$nfi <- NA
      fitMeasures$tli <- NA
      fitMeasures$nnfi <- NA
      fitMeasures$rfi <- NA
      fitMeasures$ifi <- NA
      fitMeasures$rni <- NA
      fitMeasures$cfi <- NA
  }
  
  # If LLs are not good, break here:

  if (!is.finite(Tm) ){
    return(x)
  }
  
  # RMSEA

  # fitMeasures$rmsea <- sqrt( max(Tm - dfm,0) / (sampleSize * dfm))
  fitMeasures$rmsea <-  sqrt( max( c((Tm/sampleSize)/dfm - 1/sampleSize, 0) ) )
  if (!is.finite(fitMeasures$rmsea)) fitMeasures$rmsea <- NA
  
  # FIXME: Multi-group correction from lavaan source code:
  nGroups <- nrow(x@sample@groups)
  fitMeasures$rmsea <-  fitMeasures$rmsea  * sqrt(nGroups)

  # Codes for rmsea confidence interval taken from lavaan:
  lower.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.95)
  }
  if(is.na(Tm) || is.na(dfm)) {
    fitMeasures$rmsea.ci.lower <- NA
  } else if(dfm < 1 || lower.lambda(0) < 0.0) {
    fitMeasures$rmsea.ci.lower <- 0
  } else {
    if (lower.lambda(0) * lower.lambda(Tm) > 0){
      lambda.l <- NA
    } else {
      lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=Tm)$root,
                      silent=TRUE)      
    }
    fitMeasures$rmsea.ci.lower <- sqrt( lambda.l/(sampleSize*dfm) ) * sqrt(nGroups)
  }
  
  N.RMSEA <- max(sampleSize, Tm*4) 
  upper.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.05)
  }
  if(is.na(Tm) || is.na(dfm)) {
    fitMeasures$rmsea.ci.upper <- NA
  } else if(dfm < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
    fitMeasures$rmsea.ci.upper <- 0
  } else {
    
    if (upper.lambda(0) * upper.lambda(N.RMSEA) > 0){
      lambda.u <- NA
    } else {
      
      lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                      silent=TRUE)  
    }
    
    if(inherits(lambda.u, "try-error")) {lambda.u <- NA }
    
    fitMeasures$rmsea.ci.upper <- sqrt( lambda.u/(sampleSize*dfm) )  * sqrt(nGroups)
  }
  
  fitMeasures$rmsea.pvalue <- 
    1 - pchisq(Tm, df=dfm, ncp=(sampleSize*dfm*0.05^2/nGroups))
  
  # information criteria:


  # Deviance based AIC (traditional definition)
  fitMeasures$aic.ll <-  -2*LL + 2* fitMeasures$npar
  # Deviance based AIC with sample size adjustment
  fitMeasures$aic.ll2 <-  -2*LL + 2* fitMeasures$npar +
    (2*fitMeasures$npar^2 + 2*fitMeasures$npar)/(sampleSize - fitMeasures$npar - 1)

  # Chi-square based AIC with df penalty (Kaplan, 2000): AIC(null) - AIC(saturated)
  fitMeasures$aic.x <- Tm - 2*fitMeasures$df

  # Chi-square based AIC with parameter penalty (Kline, 2016) - couldn't find the proper derivation
  fitMeasures$aic.x2 <- Tm + 2*fitMeasures$npar

  BIC <- -2*LL + fitMeasures$npar * log(sampleSize)
  fitMeasures$bic <- BIC

  # add sample-size adjusted bic
  N.star <- (sampleSize + 2) / 24
  BIC2 <- -2*LL + fitMeasures$npar * log(N.star)
  fitMeasures$bic2 <- BIC2

  # Add extended bic:
  fitMeasures$ebic.25 <-  -2*LL + fitMeasures$npar * log(sampleSize) + 4 *  fitMeasures$npar * 0.25 * log(nVar)
  fitMeasures$ebic.5 <-  -2*LL + fitMeasures$npar * log(sampleSize) + 4 *  fitMeasures$npar * 0.5 * log(nVar)
  fitMeasures$ebic.75 <-  -2*LL + fitMeasures$npar * log(sampleSize) + 4 *  fitMeasures$npar * 0.7 * log(nVar)
  fitMeasures$ebic1 <-  -2*LL + fitMeasures$npar * log(sampleSize) + 4 *  fitMeasures$npar * 1 * log(nVar)
  # fitMeasures$ebicTuning <- ebicTuning

  # Put in objet:
  x@fitmeasures <- fitMeasures
  return(x)
}
# 
# print.precisionFit <- function(x,...){
#   name <- deparse(substitute(x))[[1]]
#   if (nchar(name) > 10) name <- "object"
#   if (name=="x") name <- "object"
#   
#   cat("\nprecisionFit object:\n",
#       paste0("Use plot(",name,") to plot the network structure"),
#       "\n",
#       paste0("Fit measures stored under ",name,"$fitMeasures"),
#       "\n\n"
#   )
#   
#   fit <- data.frame(Measure = names(x$fitMeasures),
#                     Value = goodNum(unlist(x$fitMeasures)))
#   rownames(fit) <- NULL
#   print(fit)
# }
# 
# plot.precisionFit <- function(x,...){
#   qgraph::qgraph(x$network,...)
# }
