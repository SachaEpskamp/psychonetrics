# General-purpose parallel loop for psychonetrics
# Replaces boilerplate: makeCluster / clusterExport / parLapply / stopCluster
loop_psychonetrics <- function(
    expr,
    reps = 1000,
    nCores = 1,
    export = character(0),
    packages = c("psychonetrics", "dplyr"),
    verbose = TRUE,
    envir = parent.frame()
){
  # Capture the expression:
  expr_sub <- substitute(expr)

  # Auto-detect variables referenced in the expression that exist in envir:
  referenced <- all.vars(expr_sub)
  if (length(referenced) > 0) {
    keep <- vapply(referenced, function(v) {
      exists(v, envir = envir, inherits = TRUE) &&
        !is.function(get(v, envir = envir))
    }, logical(1))
    auto_export <- referenced[keep]
  } else {
    auto_export <- character(0)
  }
  all_exports <- unique(c(export, auto_export))

  if (nCores > 1) {
    # Create cluster:
    cl <- makePSOCKcluster(nCores)
    on.exit(stopCluster(cl), add = TRUE)

    # Load packages on workers:
    clusterCall(cl, function(pkgs) {
      for (pkg in pkgs) {
        suppressMessages(library(pkg, character.only = TRUE))
      }
    }, pkgs = packages)

    # Export variables to workers:
    if (length(all_exports) > 0) {
      clusterExport(cl, all_exports, envir = envir)
    }

    # Run with progress bar:
    results <- pblapply(seq_len(reps), function(i) {
      eval(expr_sub)
    }, cl = cl)
  } else {
    # Sequential with progress bar:
    results <- pblapply(seq_len(reps), function(i) {
      eval(expr_sub, envir = envir)
    })
  }

  return(results)
}
