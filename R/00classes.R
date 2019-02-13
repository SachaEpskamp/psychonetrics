# psychonetrics sample:
generate_psychonetrics_samplestats <- setClass("psychonetrics_samplestats",  slots = c(
  covmat = "array",
  means = "matrix",
  groups = "data.frame", # Data frame with information on each group  
  variables = "data.frame"
), prototype = list(groups = data.frame(
  label = character(0),
  id = integer(0),
  nobs = integer(0),stringsAsFactors = FALSE
),
variables = data.frame(
  label = character(0),
  id = integer(0)
)))

# Parameter table:


# Psychonetrics model:
generate_psychonetrics <- setClass("psychonetrics", slots = c(
  model = "character", # Model framework
  parameters = "data.frame", # Parameter table data.frame(from,  edge, to,  est,  std,  se,  matrix,  row,  col,  par)
  computed = "logical", # Logical, is the model computed yet?
  sample = "psychonetrics_samplestats", # Sample statistics
  extraMatrices = "list" # General extra matrices list
  ),
  prototype = list(
    model = "dummy",
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
      matrix = character(0),
      row = numeric(0),
      col = numeric(0),
      par = integer(0),
      group = character(0),
      group_id = integer(0),
      fixed = logical(0),stringsAsFactors = FALSE
    ),
    computed = FALSE
  ))

generate_psychonetrics()
