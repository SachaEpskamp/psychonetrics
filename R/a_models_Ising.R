# Latent network model creator
Ising <- function(
  data, # Dataset
  omega = "full", # Partial correlations
  tau,
  beta,
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
  optimizer = "default",
  storedata = FALSE,
  WLS.V,
  sampleStats, # Leave to missing
  identify = TRUE,
  verbose = TRUE,
  maxNodes = 20
){
  covtype <- match.arg(covtype)

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

    if (length(responses) != 2){
      stop("Only binary responses that are encoded in the same way are supported.")
    }
  }

  
  if (estimator == "default"){
    # if (length(ordered) > 0){
    #   estimator <- "WLS"
    # } else {
      estimator <- "ML"
    # }
  }
  
  # Fail if estimator is not ML (nothing else supported yet):
  if (estimator != "ML"){
    stop("Only ML estimation is currently supported for Ising model.")
  }
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    # WLS weights:
    if (missing(WLS.V)){
      WLS.V <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
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
                               weightsmatrix = WLS.V,
                               # corinput = corinput,
                               covtype=covtype,
                               verbose=verbose)
  }
 
  
 
  # Check some things:
  nNode <- nrow(sampleStats@variables)
  
  # Check number of nodes:
  if (nNode > maxNodes){
    stop("Aborting because the number of nodes is larger than 'maxNodes'. High-dimensional Ising models are not possible with ML estimation.")
  }
  
  # Make type:
  type <- paste0("(",responses[1]," & ",responses[2],")")
  
  # Generate model object:
  model <- generate_psychonetrics(model = "Ising",sample = sampleStats, computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Ising",
                                  rawts = FALSE, types = list(),
                                  submodel = type)
  
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
                                         sampletable = sampleStats)
  

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
    model@extramatrices <- list(
      D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
      L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
      In = Diagonal(nNode), # Identity of dim n
      responses = responses
    )

  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # Form baseline model:
    model@baseline_saturated$baseline <- Ising(data,
                                                  omega = "empty",
                                                  vars = vars,
                                                  groups = groups,
                                                  covs = covs,
                                                  means = means,
                                                  nobs = nobs,
                                                  missing = missing,
                                                  equal = equal,
                                                  estimator = estimator,
                                                  responses = responses,
                                                  baseline_saturated = FALSE,sampleStats=sampleStats)
    
    # Identify model:
    if (identify){
      model@baseline_saturated$baseline <- identify(model@baseline_saturated$baseline)
    }
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- Ising(data,
                                               omega = "full",
                                               vars = vars,
                                               groups = groups,
                                               covs = covs,
                                               means = means,
                                               nobs = nobs,
                                               missing = missing,
                                               equal = equal,
                                               estimator = estimator,
                                               responses = responses,
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
  
  # Return model:
  return(model)
}