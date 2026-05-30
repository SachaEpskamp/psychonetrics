# Compute Z:
# Generalized to any number of ordered response options, encoded identically
# across all variables, and to the Spin distribution's per-node quadratic
# (Blume-Capel) term delta. For length-2 'responses' and delta = 0 this returns
# the same partition function as the original two-outcome Ising implementation.
computeZ <- function(graph, thresholds, beta, responses = c(-1, 1), min_sum = -Inf, delta){
  graph <- as.matrix(graph)
  stopifnot(isSymmetric(unname(graph)))
  if (any(diag(graph) != 0)) {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  if (missing(delta)){
    delta <- rep(0, nrow(graph))
  }
  computeZ_cpp(graph, as.vector(thresholds), as.vector(delta), as.numeric(beta),
               as.numeric(responses), min_sum)
}
