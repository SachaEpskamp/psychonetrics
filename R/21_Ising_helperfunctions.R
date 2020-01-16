# Compute Z:
computeZ <- function(graph,thresholds,beta,responses = c(-1,1)){
  # library("IsingSampler")
  stopifnot(isSymmetric(graph))
  stopifnot(length(responses) == 2)
  if (any(diag(graph) != 0)) {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  N <- nrow(graph)
  Allstates <- do.call(expand.grid, lapply(1:N, function(x) c(responses[1], 
                                                              responses[2])))
  P <- exp(-beta * apply(Allstates, 1, function(s) IsingSampler:::H(graph, 
                                                                    s, thresholds)))
  sum(P)
}
