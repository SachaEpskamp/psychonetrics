# psychonetrics sample:
generate_psychonetrics_samplestats <- setClass("psychonetrics_samplestats",  slots = c(
  covs = "list",
  cors = "list",
  means = "list",
  thresholds = "list",
  squares = "list",
  groups = "data.frame", # Data frame with information on each group  
  variables = "data.frame",
  nobs = "numeric", # Number of observations
  corinput = "logical",
  # missingness = "list", # Missing patterns, only used when rawts = TRUE
  # data = "list" # Raw data, used only with fimldata
  fimldata = "list",
  fullFIML = "logical",
  WLS.W = "list", # List with weights matrix per group
  rawdata = "data.frame" # For bootstrapping!
), prototype = list(groups = data.frame(
  label = character(0),
  id = integer(0),
  nobs = integer(0),
  stringsAsFactors = FALSE
),
variables = data.frame(
  label = character(0),
  id = integer(0),
  ordered = logical(0)
),
corinput = FALSE,
fullFIML = FALSE
))

# Timestamp:
setOldClass("sessionInfo")

psychonetrics_log <- setClass("psychonetrics_log",  slots = c(
  event = "character",
  time = "POSIXct",
  sessionInfo = "sessionInfo"))


generate_psychonetrics_logentry <- function(event){
  stamp <- psychonetrics_log()
  stamp@event <- event
  stamp@time <- Sys.time()
  # stamp@sessionInfo <- sessionInfo()
  stamp
}

addLog <- function(x,event){
  x@log[[length(x@log)+1]] <- generate_psychonetrics_logentry(event)
  x
}

createLogList <- function(){
  res <- list(generate_psychonetrics_logentry("Model created"))
  class(res) <- "psychonetrics_log"
  return(res)
}

# Psychonetrics model:
generate_psychonetrics <- setClass("psychonetrics", slots = c(
  model = "character", # Model framework
  submodel = "character",
  parameters = "data.frame", # Parameter table data.frame(from,  edge, to,  est,  std,  se,  matrix,  row,  col,  par)
  matrices = "data.frame",
  computed = "logical", # Logical, is the model computed yet?
  sample = "psychonetrics_samplestats", # Sample statistics
  modelmatrices = "list", # Model matrices in list form
  # fitfunctions = "list", # contains fitfunction, gradient, hessian and extramatrices, logliks
  log = "psychonetrics_log",
  optim = "list",
  fitmeasures = "list",
  baseline_saturated = "list",
  equal = "character",
  objective = "numeric",
  information = "matrix",
  identification = "character",
  optimizer = "character",
  estimator = "character",
  distribution = "character",
  extramatrices = "list", # Contains extra matrices
  rawts = "logical",
  Drawts = "list",
  types = "list",
  cpp = "logical",
  meanstructure = "logical",
  verbose = "logical"
),
prototype = list(
  model = "dummy", submodel = "none",
  parameters = data.frame(
    var1 = character(0),
    var1_id = integer(0),
    op = character(0),
    var2 = character(0),
    var2_id = integer(0),
    est = numeric(0),
    std = numeric(0),
    se = numeric(0),
    p = numeric(0),
    se_boot = numeric(0),
    p_boot = numeric(0),
    matrix = character(0),
    row = numeric(0),
    col = numeric(0),
    par = integer(0),
    group = character(0),
    group_id = integer(0),
    fixed = logical(0),
    symmetrical = logical(0), # Used to determine if matrix is symmetrical
    mi = numeric(0), # Modification index
    pmi = numeric(0), #p-value modification index
    epc = numeric(0),
    mi_free = numeric(0), # Modification index
    pmi_free = numeric(0), #p-value modification index
    epc_free = numeric(0),
    mi_equal = numeric(0), # Modification index constraning groups to be equal
    pmi_equal = numeric(0), #p-value modification index constraining groups to be equal
    pmi_free = numeric(0), #p-value modification index constraining groups to be equal
    minimum = numeric(0),
    maximum = numeric(0),
    identified = logical(0), # Indicating a parameter is fixed to identify the model!
    stringsAsFactors = FALSE
  ),
  matrices = data.frame(
    name = character(0),
    nrow = integer(0),
    ncol = integer(0),
    ngroup = integer(0),
    symmetrical = logical(0),
    sparse = logical(0),
    posdef = logical(0),
    diagonal = logical(0),
    incomplete = logical(0)
  ),
  computed = FALSE,
  log = createLogList(),
  identification = "none",
  optimizer = "ucminf",
  estimator = "ML",
  rawts = FALSE,
  cpp = TRUE, # Use C++ when available
  meanstructure = TRUE,
  verbose = FALSE
))

# generate_psychonetrics()
