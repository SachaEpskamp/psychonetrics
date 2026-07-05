# Compare function for psychonetrics models:
compare <- function(...,
                    scaled.test.method = c("auto",
                                           "satorra.bentler.2001",
                                           "satorra.bentler.2010",
                                           "satorra.2000")){
  # Obtain dots:
  dots <- list(...)

  scaled.test.method <- match.arg(scaled.test.method)

  # if anything is not a psychonetrics object, stop:
  classes <- sapply(dots, class)
  if (!all(classes == "psychonetrics")){
    stop("psychonetrics::compare(...) must only contain psychonetrics objects")
  }
  
  # If not all computed, stop:
  if (!all(sapply(dots,function(x)x@computed))){
    stop("Not all models are computed.")
  }
  
  # If names are null, add names:
  if (is.null(names(dots))){
    names(dots) <- paste0("Model ",seq_along(dots))
  }
  
  # If there are not at least two models, compare to saturated:
  if (length(dots)==1){
    stop("psychonetrics::compare(...) requires at least two models")  
    # dots <- c(list(saturated = dots[[1]]@baseline_saturated$saturated),dots,list(baseline = dots[[1]]@baseline_saturated$saturated))
  }
  
  # Obtain the fitmeasures of each model:
  fits <- lapply(dots, function(x) x@fitmeasures)

  # Create first table:
  Tab <- data.frame(
    model = names(dots),
    DF = sapply(fits,"[[","df"),
    AIC = sapply(fits,"[[","aic.ll"),
    BIC = sapply(fits,"[[","bic"),
    RMSEA = sapply(fits,"[[","rmsea"),
    Chisq = sapply(fits,"[[","chisq"),
    Chisq_diff = NA,
    DF_diff = NA,
    p_value = NA
  )
  
  # Arrange table by Df:
  Tab <- Tab %>% arrange(.data[['DF']])
  
  # Compute chi-square difference (signed: should be >= 0 for nested models,
  # a negative value flags a non-nested or improperly ordered comparison):
  Tab$Chisq_diff <- c(NA,diff(Tab$Chisq))

  # Compute DF diff:
  Tab$DF_diff <- c(NA,diff(Tab$DF))

  # Compute p (only meaningful for nested models with a non-zero df difference):
  Tab$p_value <- pchisq(Tab$Chisq_diff,Tab$DF_diff,lower.tail=FALSE)

  # Models with equal df cannot be compared with a chi-square difference test:
  Tab$p_value[!is.na(Tab$DF_diff) & Tab$DF_diff == 0] <- NA

  # Warn for negative chi-square differences (non-nested or mis-ordered models):
  if (any(!is.na(Tab$Chisq_diff) & Tab$Chisq_diff < 0)){
    warning("Negative chi-square difference encountered; models may not be nested. Chi-square difference test is not valid.")
  }

  # Satorra-Bentler-family SCALED chi-square difference test. Only added when
  # ALL compared models are fit with a robust estimator (MLM/MLMV/MLMVS/MLR or
  # DWLS/WLS/ULS). For ordinary ML/FIML/etc. the scaled columns are absent so
  # that the output is byte-identical to earlier versions (backward compat).
  # A model qualifies if it carries a finite scaling factor, OR it is a saturated
  # robust model (df == 0): a saturated model has no scaling factor of its own
  # (the scaling divides by df) but can still serve as the less-constrained
  # reference, where the scaling factor is taken to be 1 (Satorra & Bentler 2001).
  is_robust_estimator <- function(x){
    if (x@estimator %in% c("WLS", "DWLS", "ULS")) return(TRUE)
    cfg <- tryCatch(get_robust_config(x), error = function(e) list())
    isTRUE(nzchar(cfg$se)) || isTRUE(nzchar(cfg$test))
  }
  has_scaling <- function(x){
    sf <- x@fitmeasures$chisq.scaling.factor
    if (!is.null(sf) && length(sf) == 1 && is.finite(sf)) return(TRUE)
    # Saturated robust model (df 0): no scaling factor, still eligible as M1.
    df_x <- x@fitmeasures$df
    is_robust_estimator(x) && !is.null(df_x) && length(df_x) == 1 && df_x == 0
  }
  if (length(dots) >= 2 && all(vapply(dots, has_scaling, logical(1)))){
    # The dots, arranged by DF ascending to match the table ordering:
    dots_sorted <- dots[order(sapply(fits, "[[", "df"))]

    # Resolve scaled.test.method = "auto" to the correct difference test for the
    # models' scaled test TYPE, exactly as lavaan's lavTestLRT does. The naive
    # Satorra-Bentler 2001 mean-scaling is only valid for the mean-adjusted tests
    # (MLM / MLR / WLSM); the scaled-shifted (MLMV / WLSMV) and mean-and-variance
    # adjusted (MLMVS) estimators require the exact Satorra (2000) trace-based
    # difference test (scaled-and-shifted, resp. mean-adjusted). Applying SB-2001
    # to a scaled-shifted fit silently returns a WRONG statistic, so "auto" is
    # the default:
    scaled_type <- function(x){
      if (x@estimator %in% c("WLS", "DWLS", "ULS")) return("scaled.shifted")
      cfg <- tryCatch(get_robust_config(x), error = function(e) list())
      tst <- if (is.null(cfg$test) || !nzchar(cfg$test)) "" else cfg$test
      switch(tst,
        "scaled.shifted"    = "scaled.shifted",
        "mean.var.adjusted" = "mean.var",
        "mean")   # satorra.bentler, yuan.bentler.mplus, or none -> mean adjusted
    }
    if (identical(scaled.test.method, "auto")){
      types <- vapply(dots_sorted, scaled_type, character(1))
      # lavaan's lavTestLRT default: the exact Satorra (2000) scaled-and-shifted
      # difference test for BOTH the scaled-shifted (MLMV/WLSMV) and the
      # mean-and-variance adjusted (MLMVS) estimators, and the Satorra-Bentler
      # 2001 mean-scaling for the mean-adjusted tests (MLM / MLR / WLSM).
      if (any(types %in% c("scaled.shifted", "mean.var"))){
        scaled.test.method <- "satorra.2000"; scaled.shifted <- TRUE
      } else {
        scaled.test.method <- "satorra.bentler.2001"; scaled.shifted <- FALSE
      }
    } else {
      # An explicitly requested Satorra (2000) test is reported in its
      # scaled-and-shifted form (the lavaan default for that method).
      scaled.shifted <- identical(scaled.test.method, "satorra.2000")
    }

    Tab$Chisq_diff_scaled <- NA_real_
    Tab$p_value_scaled <- NA_real_
    for (i in seq_len(nrow(Tab) - 1L)){
      m1 <- dots_sorted[[i]]        # less constrained (lower df)
      m0 <- dots_sorted[[i + 1L]]   # more constrained (higher df), nested in m1
      # Only a properly ordered, non-equal-df pair is testable:
      if (is.na(Tab$DF_diff[i + 1L]) || Tab$DF_diff[i + 1L] <= 0) next
      res <- scaled_diff_test_pair(m1, m0, method = scaled.test.method,
                                   scaled.shifted = scaled.shifted)
      Tab$Chisq_diff_scaled[i + 1L] <- res$stat
      Tab$p_value_scaled[i + 1L] <- res$pvalue
    }
    attr(Tab, "scaled.test.method") <- scaled.test.method
    attr(Tab, "scaled.shifted") <- scaled.shifted
  }

  # Set saturated chisquare to NA:
  if (any(Tab$model == "saturated"))
  Tab$Chisq[Tab$model == "saturated"] <- NA

  class(Tab) <- c("psychonetrics_compare","data.frame")

  return(Tab)
}

# Nice print function:
print.psychonetrics_compare <- function(x,...){

  scaled_method <- attr(x, "scaled.test.method")

  x$AIC <- goodNum2(x$AIC)
  x$BIC <- goodNum2(x$BIC)
  x$RMSEA <- goodNum2(x$RMSEA)
  x$Chisq <- goodNum2(x$Chisq)
  x$Chisq_diff <- goodNum2(x$Chisq_diff)
  x$DF_diff <- goodNumInt(x$DF_diff)
  x$p_value <- goodNum(x$p_value)
  if (!is.null(x$Chisq_diff_scaled)) x$Chisq_diff_scaled <- goodNum2(x$Chisq_diff_scaled)
  if (!is.null(x$p_value_scaled))    x$p_value_scaled    <- goodNum(x$p_value_scaled)


  # Make all numbers nicer
  # for (i in 5:ncol(x)){
    # x[,i] <- goodNum(x[,i])
  # }

  # Output something:
  # cat(
  #   paste0("\t\t####################################\n",
  #          "\t\t## psychonetrics model comparison ##\n",
  #          "\t\t####################################\n\n"))
  # psychonetrics_print_logo()
  # No awesome header :(

  print.data.frame(x,row.names=FALSE)

  cat("\nNote: Chi-square difference test assumes models are nested.")

  if (!is.null(scaled_method)){
    scaled_shifted <- attr(x, "scaled.shifted")
    method_label <- switch(scaled_method,
      "satorra.bentler.2001" = "Satorra & Bentler (2001)",
      "satorra.bentler.2010" = "Satorra & Bentler (2010)",
      "satorra.2000"         = if (isTRUE(scaled_shifted)) "Satorra (2000), scaled and shifted"
                               else "Satorra (2000), mean adjusted",
      scaled_method)
    cat("\nNote: Chisq_diff_scaled / p_value_scaled use the ",
        method_label, " scaled chi-square difference test.", sep = "")
  }
}
