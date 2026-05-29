# Compute Z:
# Generalized to any number of ordered (integer) response options, encoded
# identically across all variables. For length-2 'responses' this returns the
# same partition function as the original two-outcome implementation.
computeZ <- function(graph, thresholds, beta, responses = c(-1, 1), min_sum = -Inf){
  graph <- as.matrix(graph)
  stopifnot(isSymmetric(unname(graph)))
  if (any(diag(graph) != 0)) {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  computeZ_cpp(graph, as.vector(thresholds), as.numeric(beta),
               as.numeric(responses), min_sum)
}
