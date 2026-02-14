# Shared helper to standardize input arguments across model families.
# Called early in each constructor with NULL-converted arguments.
standardize_input <- function(
  data = NULL,
  covs = NULL,
  cors = NULL,
  nobs = NULL,
  corinput = NULL,
  groups = NULL,
  groupvar = NULL,
  studyvar = NULL,
  vars = NULL,
  family = "varcov",
  is_meta = FALSE,
  estimator = "ML",
  caller = "model()"
){

  # --- 1. groupvar / groups resolution ---
  if (!is.null(groupvar) && !is.null(groups)){
    warning("Both 'groupvar' and 'groups' were supplied. 'groupvar' takes precedence. ",
            "'groups' is deprecated; please use 'groupvar' instead.",
            call. = FALSE)
    groups <- groupvar
  } else if (!is.null(groupvar)){
    groups <- groupvar
  }
  # If only groups is supplied, keep it as-is (backward compat, no warning)

  # Meta families do not support multi-group:
  if (is_meta && !is.null(groups)){
    stop("Multi-group support is not yet included for meta-analytic models (",
         caller, ").", call. = FALSE)
  }

  # --- 2. Mutual exclusion: data vs covs/cors/nobs ---
  has_data <- !is.null(data)
  has_covs <- !is.null(covs)
  has_cors <- !is.null(cors)
  has_nobs <- !is.null(nobs)

  if (has_data && (has_covs || has_cors || has_nobs)){
    stop("'data' cannot be used together with 'covs', 'cors', or 'nobs' in ",
         caller, ". Please supply either raw data or summary statistics, not both.",
         call. = FALSE)
  }

  # --- 3. cors handling (when cors supplied, no data) ---
  if (has_cors && !has_data){
    if (!has_nobs){
      stop("'nobs' is required when using 'cors' as input in ", caller, ".",
           call. = FALSE)
    }
    if (family %in% c("varcov", "meta_varcov")){
      # For varcov/meta_varcov: default corinput to TRUE
      if (is.null(corinput)) corinput <- TRUE
      covs <- cors
    } else {
      # For all other families: treat cors as covs with a warning
      warning("Correlation matrices supplied via 'cors' are treated as covariance ",
              "matrices in ", caller, ". This is appropriate for standardized data ",
              "but may give incorrect results otherwise.",
              call. = FALSE)
      covs <- cors
      if (is.null(corinput)) corinput <- FALSE
    }
    has_covs <- TRUE
  }

  # --- 4. corinput = TRUE validation ---
  if (!is.null(corinput) && isTRUE(corinput)){
    if (!family %in% c("varcov", "meta_varcov")){
      stop("corinput = TRUE is not yet supported for ", caller,
           ". Please use covariance matrices as input.", call. = FALSE)
    }
  }

  # --- 5. data + corinput = TRUE: error for FIML ---
  if (has_data && !is.null(corinput) && isTRUE(corinput)){
    if (toupper(estimator) %in% c("FIML", "PFIML")){
      stop("corinput = TRUE is not supported with full-information estimators (FIML) in ",
           caller, ".", call. = FALSE)
    }
  }

  # --- 6. Meta studyvar handling ---
  if (is_meta && has_data){
    if (is.null(studyvar)){
      stop("'studyvar' is required when supplying raw 'data' to ", caller, ".",
           call. = FALSE)
    }
    if (!studyvar %in% names(data)){
      stop("'studyvar' column '", studyvar, "' not found in data.", call. = FALSE)
    }

    # Determine which columns to use for computing matrices:
    drop_cols <- studyvar
    if (!is.null(vars)){
      keep_cols <- vars
    } else {
      keep_cols <- setdiff(names(data), drop_cols)
    }

    studies <- unique(data[[studyvar]])
    nobs_out <- numeric(length(studies))
    covs_out <- vector("list", length(studies))

    use_cor <- !is.null(corinput) && isTRUE(corinput)
    if (is.null(corinput)){
      # Default for meta families:
      use_cor <- (family == "meta_varcov")
      corinput <- use_cor
    }

    for (i in seq_along(studies)){
      sub <- data[data[[studyvar]] == studies[i], keep_cols, drop = FALSE]
      nobs_out[i] <- nrow(sub)
      if (use_cor){
        covs_out[[i]] <- cor(sub, use = "pairwise.complete.obs")
      } else {
        covs_out[[i]] <- cov(sub, use = "pairwise.complete.obs")
      }
    }

    covs <- covs_out
    nobs <- nobs_out
    data <- NULL
  }

  # --- Return resolved values ---
  list(
    data = data,
    covs = covs,
    nobs = nobs,
    corinput = corinput,
    groups = groups
  )
}
