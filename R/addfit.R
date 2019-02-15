# Add fit measures to psychonetrics object!

# Nice umber function # FIXME: Not used now
goodNum <- function(x){
  sapply(x,function(xx){
    if (xx < 0.0001 & xx > 0){
      return("< 0.0001")
    }
    digits <- max(0,floor(log10(abs(xx))) + 1)
    isInt <- xx%%1 == 0
    gsub("\\.$","",formatC(signif(unlist(xx),digits=digits+(!isInt)*2), digits=digits+(!isInt)*2,format="fg", flag="#"))
  })  
}

# Computes fit measures
addfit <- function(
 x, ebicTuning = 0.25
){
  # If not computed, stop:
  if (!x@computed){
    stop("Model has not yet been run. Use runmodel(object) first!")
  }
  
  # Sample size:
  sampleSize <- sum(x@sample@groups$nobs)
  
  # FIXME: Not sure if needed still...
  mimic <- "lavaan"

  # Number of observations (not sure if this is needed, check):
  if (mimic == "lavaan"){
    Ncons <- sampleSize
  } else {
    Ncons <- sampleSize - 1
  }
 
  # Fitmeasures list:
  fitMeasures <- list()
  
  # Number of variables:
  fitMeasures$nvar <- nVar <- nrow(x@sample@variables)
  
  # Number of observations:
  fitMeasures$nobs <- x@sample@nobs
    
  # Number of parameters:
  fitMeasures$npar <- max(x@parameters$par)
  
  # Degrees of freedom:
  fitMeasures$df <- fitMeasures$nobs - fitMeasures$npar

  # Compute Fmin:
  fitMeasures$fmin <- x@objective
  fitMeasures$chisq <- 2 * Ncons * fitMeasures$fmin
  fitMeasures$pvalue <- pchisq(fitMeasures$chisq, fitMeasures$df, lower.tail = FALSE)
  
  # Baseline model:
  if (!is.null(x@baseline_saturated$baseline) && x@baseline_saturated$baseline@computed){
    fitMeasures$fmin_baseline <- x@baseline_saturated$baseline@objective
    fitMeasures$baseline.chisq <- 2 * Ncons * fitMeasures$fmin_baseline
    fitMeasures$baseline.df <- fitMeasures$nobs - max(x@baseline_saturated$baseline@parameters$par)
    fitMeasures$baseline.pvalue <- pchisq(fitMeasures$baseline.chisq, fitMeasures$baseline.df, lower.tail = FALSE)
    
    # Incremental Fit Indices
    Tb <- fitMeasures$baseline.chisq
    Tm <- fitMeasures$chisq
    
    dfb <- fitMeasures$baseline.df
    dfm <- fitMeasures$df
    
    fitMeasures$nfi <- (Tb - Tm) / Tb
    fitMeasures$tli <-  (Tb/dfb - Tm/dfm) / (Tb/dfb - 1) 
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
      fitMeasures$rfi <- NA
      fitMeasures$ifi <- NA
      fitMeasures$rni <- NA
      fitMeasures$cfi <- NA
  }
  
  
  # RMSEA
  fitMeasures$rmsea <- sqrt( max(Tm - dfm,0) / (Ncons * dfm))
  
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
    fitMeasures$rmsea.ci.lower <- sqrt( lambda.l/(sampleSize*dfm) )
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
    
    fitMeasures$rmsea.ci.upper <- sqrt( lambda.u/(sampleSize*dfm) )
  }
  
  fitMeasures$rmsea.pvalue <- 
    1 - pchisq(Tm, df=dfm, ncp=(sampleSize*dfm*0.05^2))
  
  # information criteria:

  # Saturated log-likelihood:
  satLL <- x@fitfunctions$loglik(x@baseline_saturated$saturated)
  # satLL <- ( -c -(sampleSize/2) * log(det(covMat)) - (sampleSize/2)*nVar )

  # log likelihood:
  LL <-  x@fitfunctions$loglik(x)
  
  
  fitMeasures$logl <- LL
  fitMeasures$unrestricted.logl <- satLL

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
  fitMeasures$ebic <-  -2*LL + fitMeasures$npar * log(sampleSize) + 4 *  fitMeasures$npar * ebicTuning * log(nVar)
  fitMeasures$ebicTuning <- ebicTuning

  # Put in objet:
  x@fitmeasures <- fitMeasures
  return(x)
}
# 
# print.ggmFit <- function(x,...){
#   name <- deparse(substitute(x))[[1]]
#   if (nchar(name) > 10) name <- "object"
#   if (name=="x") name <- "object"
#   
#   cat("\nggmFit object:\n",
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
# plot.ggmFit <- function(x,...){
#   qgraph::qgraph(x$network,...)
# }
