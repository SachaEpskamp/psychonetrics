# Panel data lag-1 vector autoregression (panel VAR) model creator.
#
# This is the observed-variable special case of dlvm1 (lambda = I, zero
# residual variances, no latent means), implemented as its own framework:
# the implied structures (25_panelvar_implied.R) and derivatives
# (25_panelvar_derivatives.R) skip all factor-loading algebra, which makes
# estimation faster than routing through dlvm1 with dummy matrices. The
# model matrices keep their dlvm1 names (beta, sigma_zeta_within,
# sigma_zeta_between, and their typed variants) for backward compatibility;
# the mean structure is a single 'mu' matrix of observed stationary means.
panelvar <- function(
  data, # Dataset
  vars, # Design matrix (wide) or character vector of variable names (long)

  # Type:
  within_latent = c("cov","chol","prec","ggm","cor"),
  between_latent = c("cov","chol","prec","ggm","cor"),

  # Data format:
  datatype = c("auto", "wide", "long"), # Data format: auto-detects from vars

  # Long-format support:
  idvar, # Subject ID variable (required for long format)
  beepvar, # Time point / measurement occasion variable (optional for long format)

  # Standardization:
  standardize = c("none", "z", "quantile", "z_per_wave"),

  # Temporal effects:
  beta = "full",

  # Contemporaneous (within-person) effects:
  omega_zeta_within = "full",
  delta_zeta_within = "diag",
  kappa_zeta_within = "full",
  sigma_zeta_within = "full",
  lowertri_zeta_within = "full",

  # Between-person effects:
  omega_zeta_between = "full",
  delta_zeta_between = "diag",
  kappa_zeta_between = "full",
  sigma_zeta_between = "full",
  lowertri_zeta_between = "full",

  # Mean structure:
  mu,

  # The rest:
  groups, # deprecated, use groupvar instead
  groupvar, # grouping variable (character column name or vector of group names)
  covs, # alternative covs (array nvar * nvar * ngroup)
  cors, # alternative cors (not supported in panelvar; errors if used)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  corinput, # not supported in panelvar (errors if TRUE)
  start = "version2", # <- start values, "version2" (default), "version1", "simple" or a psychonetrics model
  covtype = c("choose","ML","UB"),
  missing = "auto",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "ML",
  optimizer,
  storedata = FALSE,
  verbose = FALSE,
  sampleStats,

  baseline = c("stationary_random_intercept","stationary","independence","none"),

  bootstrap = FALSE,
  boot_sub,
  boot_resample,

  # Aliases for within_latent and between_latent:
  within,
  between,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices,  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
  # Correlation parameterization ("cor" types) matrices:
  rho_zeta_within = "full",
  SD_zeta_within = "diag",
  rho_zeta_between = "full",
  SD_zeta_between = "diag",
  # Temporal parameterization: "raw" models beta directly; "PDC" models the
  # partial directed correlations (from = row, to = column) directly:
  temporal = c("raw","PDC"),
  PDC = "full" # Only used when temporal = "PDC"
){

  covtype <- match.arg(covtype)
  datatype <- match.arg(datatype)
  standardize <- match.arg(standardize)

  # Standardize input arguments:
  .data_missing <- missing(data)
  .groups_missing <- missing(groups)
  si <- standardize_input(
    data = if(.data_missing) NULL else data,
    covs = if(missing(covs)) NULL else covs,
    cors = if(missing(cors)) NULL else cors,
    nobs = if(missing(nobs)) NULL else nobs,
    corinput = if(missing(corinput)) NULL else corinput,
    groups = if(.groups_missing) NULL else groups,
    groupvar = if(missing(groupvar)) NULL else groupvar,
    family = "panelvar", caller = "panelvar()", estimator = estimator
  )
  # Only overwrite when standardize_input actually resolved a value,
  # so that the original missing() state is preserved for downstream code:
  if (!is.null(si$data)) data <- si$data
  if (!is.null(si$covs)) covs <- si$covs
  if (!is.null(si$nobs)) nobs <- si$nobs
  if (!is.null(si$groups)) groups <- si$groups

  # CRAN Check workarounds for long-format conversion:
  variable <- NULL
  value <- NULL

  # Check start:
  if (is.character(start)){
    start <- start[1]
    if (!start %in% c("version2","version1","simple")){
      stop("'start' can only be 'version1', 'version2' (default), 'simple', or a psychonetrics object")
    }
  } else if (is(start,"psychonetrics")){
    start_mod <- start
    start <- "psychonetrics"
  } else {
    stop("'start' can only be 'version1', 'version2' (default), 'simple', or a psychonetrics object")
  }

  # Handle 'within' and 'between' aliases for within_latent and between_latent:
  if (!missing(within)){
    if (!missing(within_latent) && !identical(within_latent, eval(formals(panelvar)$within_latent))){
      warning("Both 'within' and 'within_latent' were specified; using 'within_latent'.")
    } else {
      within_latent <- within
    }
  }
  if (!missing(between)){
    if (!missing(between_latent) && !identical(between_latent, eval(formals(panelvar)$between_latent))){
      warning("Both 'between' and 'between_latent' were specified; using 'between_latent'.")
    } else {
      between_latent <- between
    }
  }

  # Match args:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  temporal <- match.arg(temporal)

  # Backward compatibility: models formed via the old dlvm1-based panelvar
  # used 'nu' for the observed intercepts. Treat an equality constraint on
  # 'nu' as one on 'mu':
  if ("nu" %in% equal){
    equal <- unique(c(equal[equal != "nu"], "mu"))
  }

  # --- Auto-detect data format ---
  if (missing(vars)){
    # If vars is missing but idvar is provided, assume long format:
    if (!missing(idvar)){
      datatype <- "long"
      # Extract vars from data columns:
      if (is.matrix(data)) data <- as.data.frame(data)
      if (!is.data.frame(data)) stop("'data' must be a data frame")
      vars <- names(data)
      vars <- vars[vars != idvar]
      if (!missing(beepvar)) vars <- vars[vars != beepvar]
      if (!.groups_missing || !is.null(si$groups)) vars <- vars[vars != groups]
    } else {
      stop("'vars' argument may not be missing")
    }
  }

  if (datatype == "auto"){
    if (is.matrix(vars)){
      datatype <- "wide"
    } else if (is.character(vars) && !is.matrix(vars)){
      datatype <- "long"
    } else {
      stop("Could not auto-detect data type from 'vars'. Please specify datatype = 'wide' or 'long'.")
    }
  }

  # --- Long-to-wide conversion ---
  if (datatype == "long"){

    # Validate idvar:
    if (missing(idvar)){
      stop("'idvar' is required when datatype = 'long' or when 'vars' is a character vector.")
    }
    if (is.matrix(data)) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("'data' must be a data frame")
    if (!is.character(vars) || is.matrix(vars)){
      stop("When datatype = 'long', 'vars' must be a character vector of variable names.")
    }
    if (!idvar %in% names(data)){
      stop("'idvar' argument does not correspond to column name of 'data'")
    }

    # Remove rows with NA cluster:
    if (any(is.na(data[[idvar]]))){
      warning("Rows with NA in idvar removed.")
      data <- data[!is.na(data[[idvar]]),]
    }

    # Handle beepvar:
    if (missing(beepvar)){
      data <- data[order(data[[idvar]]),]
      beepvar <- "BEEPVAR"
      data[[beepvar]] <- unlist(tapply(data[[idvar]], data[[idvar]], seq_along))
    }

    # Validate beepvar:
    if (any(data[[beepvar]][!is.na(data[[beepvar]])] %% 1 != 0)){
      stop("'beepvar' does not encode integer values.")
    }
    if (any(tapply(data[[beepvar]], data[[idvar]], function(x) any(duplicated(x)), simplify = TRUE))){
      stop("'beepvar' contains duplicate values for one or more cases.")
    }

    # Handle groups for long data:
    if (.groups_missing && is.null(si$groups)){
      groups <- "GROUPID"
      data[[groups]] <- "fullsample"
    }

    maxInCluster <- max(data[[beepvar]])

    # Standardize in long format (before reshape):
    if (standardize == "z"){
      for (v in seq_along(vars)){
        data[[vars[v]]] <- (data[[vars[v]]] - mean(data[[vars[v]]], na.rm = TRUE)) / sd(data[[vars[v]]], na.rm = TRUE)
      }
    } else if (standardize == "quantile"){
      for (v in seq_along(vars)){
        data[[vars[v]]] <- quantiletransform(data[[vars[v]]])
      }
    }
    # Note: "z_per_wave" and "none" do nothing here; z_per_wave handled after reshape via samplestats

    # Reshape long to wide:
    datalong <- tidyr::gather(data, variable, value, vars)
    datawide <- tidyr::pivot_wider(datalong, id_cols = c(idvar, groups),
                                    values_from = "value", names_from = c("variable", beepvar))

    # Construct design matrix:
    rowVars <- vars
    colVars <- as.character(seq(maxInCluster))
    design <- matrix(NA, length(rowVars), length(colVars))
    for (i in seq_along(rowVars)){
      for (j in seq_along(colVars)){
        # Exact matching (variable names may contain regex metacharacters):
        varName <- paste0(rowVars[i], "_", colVars[j])
        whichVar <- which(names(datawide) == varName)
        if (length(whichVar) == 1){
          design[i, j] <- varName
        }
      }
    }

    # Overwrite data and vars for downstream processing:
    data <- datawide
    vars <- design

    # Prevent double-standardization (already applied in long format):
    if (standardize %in% c("z", "quantile")){
      standardize <- "none"
    }
  }

  # --- Validate vars is now a design matrix ---
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix (with rows indicating variables and columns indicating measurements) or a character vector of variable names when using long-format data.")
  }

  # List all variables to use, in order:
  allVars <- na.omit(as.vector(vars))

  # --- Wide-format standardization (global per variable across waves) ---
  if (standardize %in% c("z", "quantile") && !.data_missing && !is.null(data)){
    if (is.matrix(data)) data <- as.data.frame(data)
    for (v in seq_len(nrow(vars))){
      varCols <- na.omit(vars[v, ])
      if (length(varCols) == 0) next

      if (standardize == "z"){
        allValues <- unlist(data[, varCols, drop = FALSE])
        globalMean <- mean(allValues, na.rm = TRUE)
        globalSD <- sd(allValues, na.rm = TRUE)
        for (col in varCols){
          data[, col] <- (data[, col] - globalMean) / globalSD
        }
      } else if (standardize == "quantile"){
        # Pool all values for this variable across waves, transform, distribute back:
        allValues <- unlist(data[, varCols, drop = FALSE])
        allTransformed <- quantiletransform(allValues)
        n <- nrow(data)
        for (ci in seq_along(varCols)){
          data[, varCols[ci]] <- allTransformed[((ci - 1) * n + 1):(ci * n)]
        }
      }
    }
    # Mark as done so samplestats doesn't re-standardize:
    standardize <- "none"
  }

  # Auto-detect missing data handling:
  if (missing == "auto") {
    if (!.data_missing && !is.null(data)){
      has_missing <- any(is.na(data[, allVars, drop = FALSE]))
      if (has_missing) {
        if (estimator == "ML") {
          estimator <- "FIML"
        } else if (estimator == "PML") {
          estimator <- "PFIML"
        } else {
          # LS variants: default to listwise
          missing <- "listwise"
        }
      } else {
        missing <- "listwise"
      }
    } else {
      missing <- "listwise"
    }
  }

  # Obtain sample stats:
  if (missing(sampleStats)){
    if (.data_missing || is.null(data)){
      sampleStats <- samplestats(vars = allVars,
                                 groups = groups,
                                 covs = covs,
                                 means = means,
                                 nobs = nobs,
                                 missing  = ifelse(estimator %in% c("FIML", "PFIML"),"pairwise",missing),
                                 fimldata = estimator %in% c("FIML", "PFIML"),
                                 storedata = storedata,
                                 covtype=covtype,
                                 standardize = if (standardize == "z_per_wave") "z" else "none",
                                 verbose=verbose,
                                 weightsmatrix = ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                                                        switch(estimator,
                                                               "WLS" = "full",
                                                               "ULS" = "identity",
                                                               "DWLS" = "diag"
                                                        )),
                                 bootstrap=bootstrap,
                                 boot_sub = boot_sub,
                                 boot_resample = boot_resample)
    } else {
      sampleStats <- samplestats(data = data,
                                 vars = allVars,
                                 groups = groups,
                                 covs = covs,
                                 means = means,
                                 nobs = nobs,
                                 missing  = ifelse(estimator %in% c("FIML", "PFIML"),"pairwise",missing),
                                 fimldata = estimator %in% c("FIML", "PFIML"),
                                 storedata = storedata,
                                 covtype=covtype,
                                 standardize = if (standardize == "z_per_wave") "z" else "none",
                                 verbose=verbose,
                                 weightsmatrix = ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                                                        switch(estimator,
                                                               "WLS" = "full",
                                                               "ULS" = "identity",
                                                               "DWLS" = "diag"
                                                        )),
                                 bootstrap=bootstrap,
                                 boot_sub = boot_sub,
                                 boot_resample = boot_resample)
    }
  }

  # Design matrix:
  design <- as.matrix(1*(!is.na(vars)))

  # time per var:
  timePerVar <- as.vector(design * row(design))
  timePerVar <- timePerVar[timePerVar!=0]

  # Number of variables:
  nVar <- nrow(vars)

  # Number of measurements:
  nTime <- ncol(vars)

  # The panelvar model requires at least two measurements per person (the
  # derivatives and implied structures use lag-1 blocks, which degenerate
  # at a single wave):
  if (nTime < 2){
    stop("At least two waves/measurements per person are required: the 'vars' design matrix has only one column.")
  }

  # row names:
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("V",seq_len(nrow(vars)))
  }
  varnames <- rownames(vars)

  # col names:
  if (is.null(colnames(vars))){
    colnames(vars) <- paste0("T",seq_len(ncol(vars)))
  }
  timenames <- colnames(vars)

  # Generate model object:
  model <- generate_psychonetrics(model = "panelvar",
                                  submodel = switch(within_latent,
                                                    "ggm" = "panelgvar",
                                                    "prec" = "panelgvar",
                                                    "panelvar"
                                  ),
                                  types = list(
                                    within_latent = within_latent, between_latent = between_latent,
                                    temporal = temporal
                                  ),
                                  sample = sampleStats,computed = FALSE,
                                  equal = equal,
                                  optimizer = defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  verbose = verbose)

  # Number of groups:
  nGroup <- nrow(model@sample@groups)

  nAllVar <- length(allVars)

  # Add number of observations:
  model@sample@nobs <-
    nAllVar * (nAllVar+1) / 2 * nGroup + # Covariances per group
    nAllVar * nGroup # Means per group

  # Model matrices:
  modMatrices <- list()

  ### STARTING VALUES ###

  # Means the same for each variant:

  # Expected means (average of the observed means per variable across waves):
  expMeans <- lapply(model@sample@means, function(x)tapply(x,timePerVar,mean,na.rm=TRUE))

  # Setup mu:
  modMatrices$mu <- matrixsetup_mu(mu,nNode = nVar, nGroup = nGroup, labels = varnames,equal = "mu" %in% equal,
                                   expmeans = expMeans, sampletable = sampleStats, name = "mu")

  if (start == "psychonetrics"){

    # Makelist:
    makelist <- function(x){
      if (!is.list(x)) {
        x <- list(x)
      }
      x
    }

    # Setup within-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_within,
                                         lowertri = lowertri_zeta_within,
                                         omega = omega_zeta_within,
                                         delta = delta_zeta_within,
                                         kappa = kappa_zeta_within, rho = rho_zeta_within, SD = SD_zeta_within,
                                         type = within_latent,
                                         name= "zeta_within",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = makelist(getmatrix(start_mod, "sigma_zeta_within")),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

    # Setup Beta (raw or PDC parameterization):
    if (temporal == "raw"){
      modMatrices$beta <- matrixsetup_beta(beta,
                                           name = "beta",
                                           nNode = nVar,
                                           nGroup = nGroup,
                                           labels = varnames,
                                           start = makelist(getmatrix(start_mod, "beta")),
                                           equal = "beta" %in% equal, sampletable = sampleStats)
    } else {
      modMatrices$PDC <- matrixsetup_PDC(PDC,
                                         name = "PDC",
                                         nNode = nVar,
                                         nGroup = nGroup,
                                         labels = varnames,
                                         equal = "PDC" %in% equal, sampletable = sampleStats,
                                         betastart = makelist(getmatrix(start_mod, "beta")),
                                         expcov = makelist(getmatrix(start_mod, "sigma_zeta_within")))
    }

    # Setup between-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between, rho = rho_zeta_between, SD = SD_zeta_between,
                                         type = between_latent,
                                         name= "zeta_between",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = makelist(getmatrix(start_mod, "sigma_zeta_between")),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

  } else if (start == "version2"){

    # The variables of the first wave are:
    firstVars <- apply(vars,1,function(x)na.omit(x)[1])

    # The variables of the last wave are:
    lastVars <- apply(vars,1,function(x)na.omit(x)[length(na.omit(x))])

    # Get some estimates per group:
    prior_estimates <- list()

    # loop per group:
    for (g in seq_len(nGroup)){
      prior_estimates[[g]] <- list()

      # First estimate the stationary distribution by averaging waves:
      prior_estimates[[g]]$stationary_estimate <- matrix(0,nrow(vars),nrow(vars))

      # And the lag-1 covariance structure:
      prior_estimates[[g]]$lag1_estimate <- matrix(0,nrow(vars),nrow(vars))

      # Loop over all pairs of variables:
      for (v in seq_len(nrow(vars))){
        for (v2 in seq_len(nrow(vars))){

          ### stationary estimate ###

          # Reset the count of pairs to zero:
          count <- 0

          # Loop over all time points:
          for (t in seq_len(ncol(vars))){

            # If both variables are included, include in the estimate:
            if (!is.na(vars[v,t]) && !is.na(vars[v2,t])){
              count <- count + 1
              prior_estimates[[g]]$stationary_estimate[v,v2] <- prior_estimates[[g]]$stationary_estimate[v,v2] + sampleStats@covs[[g]][vars[v,t],vars[v2,t]]
            }
          }

          # Average the estimates:
          prior_estimates[[g]]$stationary_estimate[v,v2] <-  prior_estimates[[g]]$stationary_estimate[v,v2] / count

          ### lag-1 estimate ###

          # Reset the count of pairs to zero:
          count <- 0

          # Loop over all time points:
          for (t in seq_len(ncol(vars)-1)){

            # If both variables are included, include in the estimate:
            if (!is.na(vars[v,t+1]) && !is.na(vars[v2,t])){
              count <- count + 1
              prior_estimates[[g]]$lag1_estimate[v,v2] <- prior_estimates[[g]]$lag1_estimate[v,v2] + sampleStats@covs[[g]][vars[v,t+1],vars[v2,t]]
            }
          }

          # Average the estimates:
          prior_estimates[[g]]$lag1_estimate[v,v2] <-  prior_estimates[[g]]$lag1_estimate[v,v2] / count
        }
      }

      # Spectral shift the stationary estimate:
      prior_estimates[[g]]$stationary_estimate <- spectralshift(prior_estimates[[g]]$stationary_estimate)

      # Finally, we also estimate the largest lag difference covariance structure, which will be the closest to the between person covariance structure:
      prior_estimates[[g]]$largest_lag_estimate <- sampleStats@covs[[g]][lastVars, firstVars]

      prior_estimates[[g]]$between_cov_estimate <- spectralshift(pmin(prior_estimates[[g]]$largest_lag_estimate ,t(prior_estimates[[g]]$largest_lag_estimate )))

      # Now we take the difference as estimate for the within person lag0 covariance matrix:
      prior_estimates[[g]]$within_cov_estimate <- spectralshift(prior_estimates[[g]]$stationary_estimate -  prior_estimates[[g]]$between_cov_estimate)

      # And the lag1 estimate for within person model:
      prior_estimates[[g]]$within_lag1_estimate <- prior_estimates[[g]]$lag1_estimate -  prior_estimates[[g]]$between_cov_estimate

      # Fix NA / NaN in lag-1 estimate:
      prior_estimates[[g]]$within_lag1_estimate <- ifelse(is.finite(prior_estimates[[g]]$within_lag1_estimate ), prior_estimates[[g]]$within_lag1_estimate , 0)

      # Estimate for beta:
      prior_estimates[[g]]$beta_estimate <- prior_estimates[[g]]$within_lag1_estimate  %*% solve_symmetric(prior_estimates[[g]]$within_cov_estimate)

      # Scale by eigenvalue:
      EVs <- eigen(prior_estimates[[g]]$beta_estimate)$values
      ev_beta <- max(sqrt(Re(EVs)^2 + Im(EVs)^2))
      if (ev_beta > 0.9){
        scalar <- 0.9 / ev_beta
        prior_estimates[[g]]$beta_estimate <- scalar * prior_estimates[[g]]$beta_estimate
      }

      # Hard truncate any above 0.5:
      prior_estimates[[g]]$beta_estimate[prior_estimates[[g]]$beta_estimate > 0.5] <- 0.5
      prior_estimates[[g]]$beta_estimate[prior_estimates[[g]]$beta_estimate < -0.5] <- -0.5
    }

    # Setup within-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_within,
                                         lowertri = lowertri_zeta_within,
                                         omega = omega_zeta_within,
                                         delta = delta_zeta_within,
                                         kappa = kappa_zeta_within, rho = rho_zeta_within, SD = SD_zeta_within,
                                         type = within_latent,
                                         name= "zeta_within",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = lapply(prior_estimates,"[[","within_cov_estimate"),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

    # Setup Beta (raw or PDC parameterization):
    if (temporal == "raw"){
      modMatrices$beta <- matrixsetup_beta(beta,
                                           name = "beta",
                                           nNode = nVar,
                                           nGroup = nGroup,
                                           labels = varnames,
                                           equal = "beta" %in% equal,
                                           start = lapply(prior_estimates,"[[","beta_estimate"),
                                           sampletable = sampleStats)
    } else {
      modMatrices$PDC <- matrixsetup_PDC(PDC,
                                         name = "PDC",
                                         nNode = nVar,
                                         nGroup = nGroup,
                                         labels = varnames,
                                         equal = "PDC" %in% equal, sampletable = sampleStats,
                                         betastart = lapply(prior_estimates,"[[","beta_estimate"),
                                         expcov = lapply(prior_estimates,"[[","within_cov_estimate"))
    }

    # Setup between-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between, rho = rho_zeta_between, SD = SD_zeta_between,
                                         type = between_latent,
                                         name= "zeta_between",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = lapply(prior_estimates,"[[","between_cov_estimate"),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

  } else if (start == "version1"){
    # Starting values as implemented up to version 0.11 (with lambda = I):

    # T=1 cov structure:
    firstVars <- apply(vars,1,function(x)na.omit(x)[1])
    secondVars <- apply(vars,1,function(x)na.omit(x)[2])
    firstSigma0 <- lapply(sampleStats@covs,function(x)spectralshift(x[firstVars,firstVars]))
    firstSigma1 <- lapply(sampleStats@covs,function(x)spectralshift(x[secondVars,firstVars]))

    # If beta = 0, these sort of estimate the within and between subject covs:
    prior_bet_cov <- lapply(firstSigma1,function(x)spectralshift(0.5*(x+t(x))))
    prior_wit_cov <- lapply(seq_along(firstSigma1),function(i)spectralshift(firstSigma0[[i]] - prior_bet_cov[[i]]))

    # Setup within-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_within,lowertri = lowertri_zeta_within,omega = omega_zeta_within,delta = delta_zeta_within,kappa = kappa_zeta_within, rho = rho_zeta_within, SD = SD_zeta_within,
                                         type = within_latent,
                                         name= "zeta_within",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = prior_wit_cov,
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

    # Setup Beta (raw or PDC parameterization):
    if (temporal == "raw"){
      modMatrices$beta <- matrixsetup_beta(beta,
                                           name = "beta",
                                           nNode = nVar,
                                           nGroup = nGroup,
                                           labels = varnames,
                                           equal = "beta" %in% equal, sampletable = sampleStats)
    } else {
      modMatrices$PDC <- matrixsetup_PDC(PDC,
                                         name = "PDC",
                                         nNode = nVar,
                                         nGroup = nGroup,
                                         labels = varnames,
                                         equal = "PDC" %in% equal, sampletable = sampleStats)
    }

    # Setup between-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between, rho = rho_zeta_between, SD = SD_zeta_between,
                                         type = between_latent,
                                         name= "zeta_between",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = lapply(prior_bet_cov, function(x) x/2),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

  } else if (start == "simple"){

    # Simple starting values

    # Setup within-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma = sigma_zeta_within,
                                         lowertri = lowertri_zeta_within,
                                         omega = omega_zeta_within,
                                         delta = delta_zeta_within,
                                         kappa = kappa_zeta_within, rho = rho_zeta_within, SD = SD_zeta_within,
                                         type = within_latent,
                                         name= "zeta_within",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = lapply(1:nGroup,function(x)diag(0.5,nVar)),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

    # Setup Beta (raw or PDC parameterization):
    if (temporal == "raw"){
      modMatrices$beta <- matrixsetup_beta(beta,
                                           name = "beta",
                                           nNode = nVar,
                                           nGroup = nGroup,
                                           labels = varnames,
                                           start = lapply(1:nGroup,function(x)diag(0.1,nVar)),
                                           equal = "beta" %in% equal, sampletable = sampleStats)
    } else {
      modMatrices$PDC <- matrixsetup_PDC(PDC,
                                         name = "PDC",
                                         nNode = nVar,
                                         nGroup = nGroup,
                                         labels = varnames,
                                         equal = "PDC" %in% equal, sampletable = sampleStats,
                                         betastart = lapply(1:nGroup,function(x)diag(0.1,nVar)),
                                         expcov = lapply(1:nGroup,function(x)diag(0.5,nVar)))
    }

    # Setup between-person varcov:
    modMatrices <- c(modMatrices,
                     matrixsetup_flexcov(sigma_zeta_between,lowertri_zeta_between,omega_zeta_between,delta_zeta_between,kappa_zeta_between, rho = rho_zeta_between, SD = SD_zeta_between,
                                         type = between_latent,
                                         name= "zeta_between",
                                         sampleStats= sampleStats,
                                         equal = equal,
                                         nNode = nVar,
                                         expCov = lapply(1:nGroup,function(x)diag(0.5,nVar)),
                                         nGroup = nGroup,
                                         labels = varnames
                     ))

  }

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  model@extramatrices <- list(
    # Entire duplication matrix needed for likelihood:
    D = psychonetrics::duplicationMatrix(nAllVar),

    # Duplication matrix per variable block:
    D2 = psychonetrics::duplicationMatrix(nVar),

    # Strict duplication matrix:
    Dstar = psychonetrics::duplicationMatrix(nVar,diag = FALSE),

    # Identity matrix:
    In = as(diag(nVar),"dMatrix"),

    # Diagonalization matrix:
    A = psychonetrics::diagonalizationMatrix(nVar),

    # Commutation matrix:
    C = commutationMatrix(nVar, nVar),

    # Elimination matrix:
    L = psychonetrics::eliminationMatrix(nVar),

    design = design
  )

  # Come up with P...
  # Dummy matrices with indices:
  muDummy <- matrix(rep(1:nVar,nTime))
  sigDummy <- matrix(0,nVar,nVar)
  sigDummy[lower.tri(sigDummy,diag=TRUE)] <- max(muDummy) + seq_len(nVar*(nVar+1)/2)
  sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]

  U <- list(sigDummy)
  # Now make all lag-k blocks...
  # Form blocks:
  for (i in 1:(nTime-1)){
    U[[length(U) + 1]] <- matrix(max(unlist(U)) + seq_len(nVar^2), nVar, nVar)
  }

  allSigmas <- blockToeplitz(U)

  # Total number:
  totElements <- max(allSigmas)

  # Now subset with only observed:
  subMu <- muDummy[as.vector(design==1),,drop=FALSE]
  subSigmas <- allSigmas[as.vector(design==1),as.vector(design==1)]

  # P matrix:
  distVec <-  c(as.vector(subMu),subSigmas[lower.tri(subSigmas,diag=TRUE)])
  nTotal <- length(distVec)
  distVecrawts <- seq_along(distVec)[distVec!=0]
  distVec <- distVec[distVec!=0]
  # Now I can make the matrix:

  model@extramatrices$P <- sparseMatrix(
    i = distVecrawts, j = distVec, dims = c(nTotal, totElements)
  )

  model@extramatrices$P <- as( model@extramatrices$P, "dMatrix")

  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)

  # Baseline model:
  baseline <- match.arg(baseline)

  ### Baseline model ###
  if (is.list(baseline_saturated)){
    model@baseline_saturated <- baseline_saturated
  } else if (isTRUE(baseline_saturated)){

    # Form baseline model:
    if (baseline ==  "stationary_random_intercept"){

      # Form the baseline model:
      model@baseline_saturated$baseline <- panelvar(data,
                                                 vars = vars,
                                                 groups = groups,
                                                 covs = covs,
                                                 means = means,
                                                 nobs = nobs,
                                                 missing = missing,
                                                 equal = equal,
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats = sampleStats,

                                                 # Empty networks:
                                                 sigma_zeta_within = "diag",
                                                 sigma_zeta_between = "diag",
                                                 beta = "zero"
                                                 )

    } else if (baseline == "stationary"){

      # Form the baseline model:
      model@baseline_saturated$baseline <- panelvar(data,
                                                    vars = vars,
                                                    groups = groups,
                                                    covs = covs,
                                                    means = means,
                                                    nobs = nobs,
                                                    missing = missing,
                                                    equal = equal,
                                                    estimator = estimator,
                                                    baseline_saturated = FALSE,
                                                    sampleStats = sampleStats,

                                                    # Empty networks:
                                                    sigma_zeta_within = "diag",
                                                    sigma_zeta_between = "zero",
                                                    beta = "zero"
      )

    } else if (baseline == "independence"){

      model@baseline_saturated$baseline <- varcov(data,
                                                  type = "chol",
                                                  lowertri = "empty",
                                                  vars = allVars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  baseline_saturated = FALSE,
                                                  sampleStats = sampleStats)

    }

    ### Saturated model ###
    # No `equal = equal`: the saturated reference is unconstrained per
    # group; cross-group equality belongs to the target/baseline only.
    model@baseline_saturated$saturated <- varcov(data,
                                                 type = "chol",
                                                 lowertri = "full",
                                                 vars = allVars,
                                                 groups = groups,
                                                 covs = covs,
                                                 means = means,
                                                 nobs = nobs,
                                                 missing = missing,
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats = sampleStats)

    # if not FIML/PFIML, Treat as computed:
    if (!estimator %in% c("FIML", "PFIML")){
      model@baseline_saturated$saturated@computed <- TRUE

      # FIXME: TODO
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
    }
  }

  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # Setup PML/PFIML penalization:
  if (estimator %in% c("PML", "PFIML")) {
    model@penalty <- list(lambda = penalty_lambda, alpha = penalty_alpha)
    pen_mats <- if (missing(penalize_matrices)) defaultPenalizeMatrices(model) else penalize_matrices
    model <- penalize(model, matrix = pen_mats, lambda = penalty_lambda, log = FALSE)
    # Baseline/saturated models should use unpenalized estimator:
    base_est <- if (estimator == "PFIML") "FIML" else "ML"
    if (!is.null(model@baseline_saturated$baseline)) {
      model@baseline_saturated$baseline@estimator <- base_est
    }
    if (!is.null(model@baseline_saturated$saturated)) {
      model@baseline_saturated$saturated@estimator <- base_est
    }
  }

  # Return model:
  return(model)
}
