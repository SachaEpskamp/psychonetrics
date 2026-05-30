# Matrix setup for the Spin-distribution quadratic node potential 'delta'
# (the Blume-Capel single-ion / crystal-field term). Like the threshold vector
# 'tau' this is one parameter per node (a nNode x nGroup matrix internally),
# but it multiplies x_i^2 rather than x_i in the exponent:
#
#   P(x) propto exp( sum_i tau_i x_i - sum_i delta_i x_i^2 + sum_{i<j} omega_ij x_i x_j )
#
# A value of 0 in 'delta' fixes that element to zero (this is how the Ising
# model is obtained from the Spin distribution: all delta_i fixed at 0). A
# missing 'delta' frees every element (the default for the BlumeCapel model).
matrixsetup_spindelta <- function(
  delta, # delta argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  labels,
  equal = FALSE,
  sampletable,
  name = "delta"
){
  # Fix delta (per-node vector; 0 => fixed at zero, 1/missing => free):
  delta <- fixMu(delta, nGroup, nNode, equal)

  # For each group, form starting values (start every node at 0, i.e. Ising):
  deltaStart <- delta
  for (g in 1:nGroup){
    deltaStart[,g][] <- 0
  }

  # Form the model matrix part. delta is a real-valued node potential (it may be
  # negative -- negative delta favours the extreme categories), so no bounds:
  list(delta,
       mat = name,
       op = "~~",
       symmetrical = FALSE,
       sampletable = sampletable,
       rownames = labels,
       colnames = "1",
       start = deltaStart)
}
