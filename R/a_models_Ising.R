# Latent network model creator
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
  covtype <- match.arg(covtype)
  beta_model <- match.arg(beta_model)
  
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

  # The Ising model supports any number of ordered response options, encoded
  # identically across all variables. Only non-integer responses are an error
  # (a single unique value is degenerate and also rejected):
  if (length(responses) < 2){
    stop("At least two distinct response options are required.")
  }
  if (any(is.na(responses)) ||
      any(abs(responses - round(responses)) > sqrt(.Machine$double.eps))){
    stop("Only integer responses are supported (encoded identically across all variables).")
  }
  responses <- round(responses)
  
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
  
  # Generate model object:
  model <- generate_psychonetrics(model = "Ising", sample = sampleStats, computed = FALSE, 
                                  equal = equal,
                                  optimizer =  "nlminb", estimator = estimator, distribution = "Ising",
                                  rawts = FALSE, types = list(beta_model = beta_model),
                                  submodel = type, verbose=verbose)
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  

  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode-1) / 2 * nGroup + # Covariances per group
    nMeans
  
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
    
    # Form baseline model:
    model@baseline_saturated$baseline <- Ising(data,
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
    model@baseline_saturated$saturated <- Ising(data,
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