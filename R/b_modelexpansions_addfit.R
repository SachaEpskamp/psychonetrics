# Add fit measures to psychonetrics object!


# Compute WLSMV (mean-and-variance adjusted) scaled test statistic.
# Implements Satorra-Bentler-style correction for WLS/DWLS/ULS estimators.
# References: Satorra & Bentler (1994), Muthén (1993)
compute_wlsmv_correction <- function(x) {

  # Check that Gamma is available:
  if (length(x@sample@WLS.Gamma) == 0) {
    return(NULL)
  }

  # Prepare model:
  if (x@cpp) {
    prep <- prepareModel_cpp(parVector(x), x)
  } else {
    prep <- prepareModel(parVector(x), x)
  }

  # Get model Jacobian (block-diagonal across groups):
  if (x@cpp) {
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov_cpp,
      "lvm" = d_phi_theta_lvm_cpp,
      "var1" = d_phi_theta_var1_cpp,
      "dlvm1" = d_phi_theta_dlvm1_cpp,
      "tsdlvm1" = d_phi_theta_tsdlvm1_cpp,
      "meta_varcov" = d_phi_theta_meta_varcov_cpp,
      "Ising" = d_phi_theta_Ising_cpp,
      "ml_lvm" = d_phi_theta_ml_lvm_cpp
    )
  } else {
    modelJacobian <- switch(
      x@model,
      "varcov" = d_phi_theta_varcov,
      "lvm" = d_phi_theta_lvm,
      "var1" = d_phi_theta_var1,
      "dlvm1" = d_phi_theta_dlvm1,
      "tsdlvm1" = d_phi_theta_tsdlvm1,
      "meta_varcov" = d_phi_theta_meta_varcov,
      "Ising" = d_phi_theta_Ising,
      "ml_lvm" = d_phi_theta_ml_lvm
    )
  }

  modelPart <- modelJacobian(prep)

  # Get M matrix (maps constrained to full parameters):
  if (x@cpp) {
    M <- Mmatrix_cpp(x@parameters)
  } else {
    M <- Mmatrix(x@parameters)
  }

  # Full Delta = Jacobian %*% M (dsigma/dtheta_free)
  Delta_full <- as.matrix(modelPart %*% M)

  # Group info:
  nGroups <- nrow(x@sample@groups)
  nobs_per_group <- x@sample@groups$nobs
  nTotal <- sum(nobs_per_group)

  # Determine row ranges per group from WLS.W dimensions:
  nrows_per_group <- sapply(x@sample@WLS.W, nrow)

  # Build per-group W (scaled) and Gamma (scaled) lists, and extract Delta blocks:
  W_list <- list()
  Gamma_list <- list()
  Delta_list <- list()

  row_offset <- 0
  for (g in seq_len(nGroups)) {
    ng <- nrows_per_group[g]
    fg <- nobs_per_group[g] / nTotal  # group weight

    # Extract per-group Delta block:
    row_inds <- row_offset + seq_len(ng)
    Delta_list[[g]] <- Delta_full[row_inds, , drop = FALSE]
    row_offset <- row_offset + ng

    # Get W for this group (weight matrix used in estimation):
    W_g <- as.matrix(x@sample@WLS.W[[g]])
    if (x@estimator == "DWLS") {
      W_g <- diag(diag(W_g))
    } else if (x@estimator == "ULS") {
      W_g <- diag(nrow(W_g))
    }

    # Get Gamma for this group (full asymptotic covariance):
    Gamma_g <- as.matrix(x@sample@WLS.Gamma[[g]])

    # Scale by group weight (following lavaan convention):
    W_list[[g]] <- W_g * fg
    Gamma_list[[g]] <- Gamma_g / fg
  }

  # Accumulate DtWD across all groups FIRST (needed for multi-group):
  npar_free <- ncol(Delta_full)
  DtWD <- matrix(0, npar_free, npar_free)
  for (g in seq_len(nGroups)) {
    DtWD <- DtWD + t(Delta_list[[g]]) %*% W_list[[g]] %*% Delta_list[[g]]
  }
  E_inv <- solve(DtWD)

  # Now compute per-group UGamma and accumulate traces:
  trace_UG <- 0
  trace_UG2 <- 0

  for (g in seq_len(nGroups)) {
    # U_g = W_g - W_g Delta_g E_inv Delta_g' W_g  (projection matrix):
    WD <- W_list[[g]] %*% Delta_list[[g]]
    U_g <- W_list[[g]] - WD %*% E_inv %*% t(WD)

    # UGamma:
    UG <- U_g %*% Gamma_list[[g]]

    # Accumulate traces:
    trace_UG <- trace_UG + sum(diag(UG))
    trace_UG2 <- trace_UG2 + sum(UG * t(UG))  # = tr(UG %*% UG)
  }

  # Compute scaling:
  if (trace_UG2 <= 0 || !is.finite(trace_UG2)) {
    return(NULL)
  }

  # Get the integer df from the model:
  nobs <- x@sample@nobs
  npar <- max(x@parameters$par)
  df_integer <- nobs - npar
  if (!is.null(x@baseline_saturated$saturated)) {
    df_integer <- df_integer - (x@baseline_saturated$saturated@sample@nobs - max(x@baseline_saturated$saturated@parameters$par))
  }

  chisq_naive <- x@objective * nTotal

  # Mean-and-Variance adjusted (Satorra-Bentler style):
  df_mv <- trace_UG^2 / trace_UG2
  scaling_mv <- trace_UG / df_mv  # = trace_UG2 / trace_UG
  chisq_mv <- chisq_naive / scaling_mv

  # Scaled-shifted (Asparouhov & Muthen 2010, used by lavaan WLSMV):
  # a = sqrt(df / trace_UG2)
  # shift = df - a * trace_UG
  # T_scaled_shifted = a * T_naive + shift
  # scaling_factor = 1/a
  a <- sqrt(df_integer / trace_UG2)
  shift <- df_integer - a * trace_UG
  chisq_ss <- a * chisq_naive + shift
  scaling_ss <- 1 / a

  list(
    # Mean-variance adjusted (fractional df):
    chisq.scaled.mv = chisq_mv,
    df.scaled.mv = df_mv,
    scaling.factor.mv = scaling_mv,
    # Scaled-shifted (integer df, matches lavaan WLSMV):
    chisq.scaled = chisq_ss,
    df.scaled = df_integer,
    scaling.factor = scaling_ss,
    shift.parameter = shift,
    # Raw traces for diagnostics:
    trace.UGamma = trace_UG,
    trace.UGamma2 = trace_UG2
  )
}


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

  # WLSMV scaled test statistic (mean-and-variance adjusted):
  wlsmv_model <- NULL
  if (x@estimator %in% c("WLS","DWLS","ULS") && length(x@sample@WLS.Gamma) > 0) {
    wlsmv_model <- tryCatch(compute_wlsmv_correction(x), error = function(e) NULL)
    if (!is.null(wlsmv_model)) {
      fitMeasures$chisq.scaled <- wlsmv_model$chisq.scaled
      fitMeasures$df.scaled <- wlsmv_model$df.scaled
      fitMeasures$pvalue.scaled <- pchisq(wlsmv_model$chisq.scaled, wlsmv_model$df.scaled, lower.tail = FALSE)
      fitMeasures$chisq.scaling.factor <- wlsmv_model$scaling.factor
    }
  }
  
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
    fitMeasures$baseline.npar <- max(x@baseline_saturated$baseline@parameters$par)
    fitMeasures$baseline.df <- fitMeasures$nobs - max(x@baseline_saturated$baseline@parameters$par)
    fitMeasures$baseline.pvalue <- pchisq(fitMeasures$baseline.chisq, fitMeasures$baseline.df, lower.tail = FALSE)

    # WLSMV scaled baseline:
    wlsmv_baseline <- NULL
    if (!is.null(wlsmv_model) && x@estimator %in% c("WLS","DWLS","ULS") && length(x@sample@WLS.Gamma) > 0) {
      wlsmv_baseline <- tryCatch(compute_wlsmv_correction(x@baseline_saturated$baseline), error = function(e) NULL)
      if (!is.null(wlsmv_baseline)) {
        fitMeasures$baseline.chisq.scaled <- wlsmv_baseline$chisq.scaled
        fitMeasures$baseline.df.scaled <- wlsmv_baseline$df.scaled
        fitMeasures$baseline.pvalue.scaled <- pchisq(wlsmv_baseline$chisq.scaled, wlsmv_baseline$df.scaled, lower.tail = FALSE)
        fitMeasures$baseline.chisq.scaling.factor <- wlsmv_baseline$scaling.factor
      }
    }
    
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

    # Scaled incremental fit indices (WLSMV):
    if (!is.null(wlsmv_model) && !is.null(wlsmv_baseline)) {
      Tm_s <- fitMeasures$chisq.scaled
      dfm_s <- fitMeasures$df.scaled
      Tb_s <- fitMeasures$baseline.chisq.scaled
      dfb_s <- fitMeasures$baseline.df.scaled

      fitMeasures$nfi.scaled <- (Tb_s - Tm_s) / Tb_s
      fitMeasures$tli.scaled <- (Tb_s/dfb_s - Tm_s/dfm_s) / (Tb_s/dfb_s - 1)
      t1_s <- (Tm_s - dfm_s)*dfb_s
      t2_s <- (Tb_s - dfb_s)*dfm_s
      fitMeasures$nnfi.scaled <- ifelse(dfm_s > 0 & abs(t2_s) > 0, 1 - t1_s/t2_s, 1)
      fitMeasures$rfi.scaled <- (Tb_s/dfb_s - Tm_s/dfm_s) / (Tb_s/dfb_s)
      fitMeasures$ifi.scaled <- (Tb_s - Tm_s) / (Tb_s - dfm_s)
      fitMeasures$rni.scaled <- ((Tb_s - dfb_s) - (Tm_s - dfm_s)) / (Tb_s - dfb_s)
      fitMeasures$cfi.scaled <- ifelse(dfm_s > Tm_s, 1, 1 - (Tm_s - dfm_s)/(Tb_s - dfb_s))
    }

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

  # Scaled RMSEA (WLSMV):
  if (!is.null(wlsmv_model)) {
    Tm_s <- fitMeasures$chisq.scaled
    dfm_s <- fitMeasures$df.scaled
    fitMeasures$rmsea.scaled <- sqrt( max( c((Tm_s/sampleSize)/dfm_s - 1/sampleSize, 0) ) )
    if (!is.finite(fitMeasures$rmsea.scaled)) fitMeasures$rmsea.scaled <- NA
    fitMeasures$rmsea.scaled <- fitMeasures$rmsea.scaled * sqrt(nGroups)
  }

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
