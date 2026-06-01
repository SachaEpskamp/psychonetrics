# Shared internal engine for the Spin-distribution models (Ising and
# BlumeCapel). Not exported; called by the Ising() and BlumeCapel() wrappers
# below (and recursively for the baseline/saturated reference models).
spinModel <- function(
  data, # Dataset
  omega = "full", # Partial correlations
  tau,
  beta,
  beta_model = c("beta","log_beta"),
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  covtype = c("choose","ML","UB"),
  responses, # May not be missing if data is missing
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer,
  storedata = FALSE,
  WLS.W,
  sampleStats, # Leave to missing
  identify = TRUE,
  verbose = FALSE,
  maxNodes = 20,
  maxStates = 2^maxNodes, # Max number of response patterns enumerated by ML (length(responses)^nNode)
  min_sum = -Inf, # Used for threhsolded Ising model estimation
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices,  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
  delta, # Spin-distribution quadratic node potential. For the Ising model this is always fixed at 0.
  model_name = c("Ising","BlumeCapel") # Which Spin-distribution model is being built.
){
  covtype <- match.arg(covtype)
  beta_model <- match.arg(beta_model)
  model_name <- match.arg(model_name)

  if (missing(data) && missing(responses)){
    stop("'responses' argument may not be missing if 'data' is missing.")
  }

  
  # Determine responses:
  if (missing(responses)){
    if (missing(vars)){
      responses <- sort(unique(unlist(c(as.matrix(data)))))
    } else {
      responses <- sort(unique(unlist(c(as.matrix(data[,vars])))))
    }
  } else {
    responses <- sort(unique(responses))
  }

  # The Spin distribution supports any number of ordered response options,
  # encoded identically across all variables. The values need not be integers;
  # only a degenerate single-value response set is rejected (NAs are dropped by
  # the sort(unique()) above):
  if (length(responses) < 2){
    stop("At least two distinct response options are required.")
  }

  # Warn when the response coding departs from the conventional one for the
  # requested model. The classical Ising model is defined for binary -1/1 or
  # 0/1 codings; the Blume-Capel model for the ternary -1/0/1 coding.
  resp_equals <- function(r, target){
    length(r) == length(target) && isTRUE(all.equal(sort(r), sort(target), check.attributes = FALSE))
  }
  if (model_name == "Ising"){
    if (!(resp_equals(responses, c(-1, 1)) || resp_equals(responses, c(0, 1)))){
      warning(paste0("Responses set to '", paste(responses, collapse = ", "),
                     "'; behavior is not exactly as the classical Ising model ",
                     "(defined for binary -1/1 or 0/1 responses)."))
    }
  } else if (model_name == "BlumeCapel"){
    # The Blume-Capel quadratic term delta_i * x_i^2 is not identifiable with
    # only two response options: x_i^2 is then constant across the two states,
    # so delta is confounded with the threshold and the normalizing constant.
    if (length(responses) < 3){
      stop(paste0("The Blume-Capel model requires at least three response options; ",
                  "with only two responses (here '", paste(responses, collapse = ", "),
                  "') the quadratic 'delta' term is not identifiable. Use Ising() for binary data."))
    }
    if (!resp_equals(responses, c(-1, 0, 1))){
      warning(paste0("Responses set to '", paste(responses, collapse = ", "),
                     "'; the Blume-Capel model is conventionally defined for ",
                     "responses -1, 0, 1."))
    }
  }

  # Check minimum sum score:
  if (min_sum > -Inf){
    if (missing(vars)){
      min_sum_in_data <- min(rowSums(as.matrix(data)))
    } else {
      min_sum_in_data <- min(rowSums(as.matrix(data[,vars])))
    }
    
    if (min_sum_in_data < min_sum){
      stop("One or more sumscores in the data are lower than the threshold set using the 'min_sum' argument.")
    }
  }

  
  if (estimator == "default"){
    # if (length(ordered) > 0){
    #   estimator <- "WLS"
    # } else {
      estimator <- "ML"
    # }
  }
  
  # Fail if estimator is not ML or PML (nothing else supported yet):
  if (!estimator %in% c("ML", "PML")){
    stop("Only ML and PML estimation are currently supported for Ising model.")
  }
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    # WLS weights:
    if (missing(WLS.W)){
      WLS.W <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                      switch(estimator,
                             "WLS" = "full",
                             "ULS" = "identity",
                             "DWLS" = "diag"
                      ))
    }
  
    sampleStats <- samplestats(data = data, 
                               vars = vars, 
                               # ordered=ordered,
                               groups = groups,
                               covs = covs, 
                               means = means, 
                               nobs = nobs, 
                               missing = ifelse(estimator == "FIML","pairwise",missing),
                               fimldata = estimator == "FIML",
                               storedata = storedata,
                               weightsmatrix = WLS.W,
                               # corinput = corinput,
                               covtype=covtype,
                               verbose=verbose,
                               bootstrap=bootstrap,
                               boot_sub = boot_sub,
                               boot_resample = boot_resample)
  }
 
  # Check some things:
  nNode <- nrow(sampleStats@variables)

  # Check size of the state space. ML estimation of the Ising model enumerates
  # every possible response pattern, so the cost grows as length(responses)^nNode.
  # The default 'maxStates' = 2^maxNodes reproduces the historical binary node
  # cap (2^nNode > 2^maxNodes  <=>  nNode > maxNodes) while accounting for more
  # than two response options. Raise 'maxStates' to override.
  nStates <- length(responses)^nNode
  if (nStates > maxStates){
    stop(paste0(
      "Aborting because the number of response patterns to enumerate (",
      length(responses), "^", nNode, " = ", format(nStates, scientific = TRUE, digits = 3),
      ") exceeds 'maxStates' (", format(maxStates, scientific = TRUE, digits = 3),
      "). Exact ML estimation enumerates every response pattern, which grows as ",
      "length(responses)^nNode. Reduce the number of nodes or response options, or ",
      "raise 'maxStates' if the computation is feasible."))
  }

  # Make type:
  type <- paste0("(", paste(responses, collapse = " & "), ")")
  
  # Generate model object. Both the Ising and BlumeCapel models are built on the
  # 'Spin' distribution; they differ only in whether the quadratic 'delta'
  # parameters are free (BlumeCapel) or fixed at zero (Ising):
  model <- generate_psychonetrics(model = model_name, sample = sampleStats, computed = FALSE,
                                  equal = equal,
                                  optimizer =  "nlminb", estimator = estimator, distribution = "Spin",
                                  rawts = FALSE, types = list(beta_model = beta_model),
                                  submodel = type, verbose=verbose)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  

  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  # Add number of observed summary statistics. The Spin distribution matches
  # the means (via tau) and the pairwise products (via omega); the BlumeCapel
  # model additionally matches the per-node second moments E(x_i^2) (via delta),
  # which contributes nNode statistics per group. (The reported degrees of
  # freedom are computed relative to the saturated model, where this count
  # cancels; it only affects the reported number of summary statistics.)
  model@sample@nobs <-
    nNode * (nNode-1) / 2 * nGroup + # Pairwise products (omega) per group
    nMeans +                         # Means (tau)
    (if (model_name == "BlumeCapel") nNode * nGroup else 0) # Second moments (delta)
  
  # Model matrices:
  modMatrices <- list()
  
  # Setup thresholds:
  modMatrices$mu <- matrixsetup_isingtau(tau, nNode = nNode,nGroup = nGroup,labels = sampleStats@variables$label,equal = "tau" %in% equal,
                       sampletable = sampleStats)
  
  
  # Setup network:
  modMatrices$omega <- matrixsetup_isingomega(omega,
                               nNode = nNode,
                               nGroup = nGroup,
                               labels = sampleStats@variables$label,
                               equal = "omega" %in% equal, sampletable = sampleStats)

  # Setup quadratic node potential (delta). For the Ising model delta is forced
  # to zero (fixed); for the BlumeCapel model it defaults to free (a missing
  # 'delta' frees every node, matching matrixsetup_spindelta/fixMu). The matrix
  # is inserted before beta so the distribution parameter order is
  # (tau, omega[lower.tri], delta, beta), matching the estimator gradient,
  # expected Hessian, and d_phi_theta of the Spin distribution.
  if (model_name == "Ising"){
    delta_spec <- 0 # all elements fixed at zero
  } else {
    delta_spec <- if (missing(delta)) "default" else delta # "default" -> all free
  }
  if (identical(delta_spec, "default")){
    modMatrices$delta <- matrixsetup_spindelta(nNode = nNode, nGroup = nGroup,
                                 labels = sampleStats@variables$label,
                                 equal = "delta" %in% equal, sampletable = sampleStats)
  } else {
    modMatrices$delta <- matrixsetup_spindelta(delta_spec, nNode = nNode, nGroup = nGroup,
                                 labels = sampleStats@variables$label,
                                 equal = "delta" %in% equal, sampletable = sampleStats)
  }

  # Setup temperature:
  modMatrices$beta <- matrixsetup_isingbeta(beta, nGroup = nGroup, equal = "beta" %in% equal,
                                         sampletable = sampleStats, log = (beta_model == "log_beta"))
  

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplication matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elimination matrix
      In = as(diag(nNode),"dMatrix"), # Identity of dim n
      responses = responses,
      min_sum = min_sum
    )

  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Form baseline model (same model type and delta treatment as the target):
    model@baseline_saturated$baseline <- spinModel(data,
                                                  beta_model = beta_model,
                                                  omega = "zero",
                                                  vars = vars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  responses = responses,
                                                  maxStates = maxStates,
                                                  delta = delta_spec,
                                                  model_name = model_name,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats,
                                               min_sum=min_sum)
    
    # Identify model:
    if (identify){
      model@baseline_saturated$baseline <- identify(model@baseline_saturated$baseline)
    }
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    # No `equal = equal`: the saturated reference is unconstrained per
    # group; cross-group equality belongs to the target/baseline only.
    model@baseline_saturated$saturated <- spinModel(data,
                                               beta_model = beta_model,
                                               omega = "full",
                                               vars = vars,
                                               groups = groups,
                                               covs = covs,
                                               means = means,
                                               nobs = nobs,
                                               missing = missing,
                                               estimator = estimator,
                                               responses = responses,
                                               maxStates = maxStates,
                                               delta = delta_spec,
                                               model_name = model_name,
                                               baseline_saturated = FALSE,sampleStats=sampleStats)
    # Identify model:
    if (identify){
      model@baseline_saturated$saturated <- identify(model@baseline_saturated$saturated)
    }
  }
  
  # Identify model:
  if (identify){
    model <- identify(model)
  }
  
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # Setup PML penalization:
  if (estimator == "PML") {
    model@penalty <- list(lambda = penalty_lambda, alpha = penalty_alpha)
    pen_mats <- if (missing(penalize_matrices)) defaultPenalizeMatrices(model) else penalize_matrices
    model <- penalize(model, matrix = pen_mats, lambda = penalty_lambda, log = FALSE)
    # Baseline/saturated models should use ML, not PML:
    if (!is.null(model@baseline_saturated$baseline)) {
      model@baseline_saturated$baseline@estimator <- "ML"
    }
    if (!is.null(model@baseline_saturated$saturated)) {
      model@baseline_saturated$saturated@estimator <- "ML"
    }
  }

  # Return model:
  return(model)
}

# Ising model creator. Thin wrapper around the shared Spin-distribution engine
# spinModel() with the quadratic 'delta' parameters fixed at zero, reproducing
# the classical Ising model exactly. See BlumeCapel() for the model that frees
# delta. Arguments supplied by the user are forwarded verbatim via match.call(),
# so the engine's own defaults and missing()-handling apply to anything omitted.
Ising <- function(
  data, # Dataset
  omega = "full", # Partial correlations
  tau,
  beta,
  beta_model = c("beta","log_beta"),
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  covtype = c("choose","ML","UB"),
  responses, # May not be missing if data is missing
  missing = "listwise",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer,
  storedata = FALSE,
  WLS.W,
  sampleStats, # Leave to missing
  identify = TRUE,
  verbose = FALSE,
  maxNodes = 20,
  maxStates = 2^maxNodes, # Max number of response patterns enumerated by ML (length(responses)^nNode)
  min_sum = -Inf, # Used for threhsolded Ising model estimation
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments:
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
){
  # Use the function object itself (not the name) in the call head: spinModel is
  # an internal, non-exported function and would not be found by name lookup
  # from the caller's environment when the package is installed.
  mc <- match.call()
  mc[[1L]] <- spinModel
  mc$model_name <- "Ising"
  eval(mc, parent.frame())
}
