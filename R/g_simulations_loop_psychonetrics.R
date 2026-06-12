# General-purpose parallel loop for psychonetrics
# Replaces boilerplate: makeCluster / clusterExport / parLapply / stopCluster
loop_psychonetrics <- function(
    expr,
    reps = 1000,
    nCores = 1,
    export = character(0),
    packages = c("psychonetrics", "dplyr"),
    verbose = TRUE,
    envir = parent.frame(),
    seed = NULL # If not NULL, used for reproducible RNG (clusterSetRNGStream in parallel; set.seed sequentially)
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

  # Also auto-export user-defined functions called in the expression. Only
  # functions defined directly in 'envir' are exported (inherits = FALSE), so
  # package/base functions are never copied to the workers:
  fun_names <- setdiff(all.names(expr_sub, unique = TRUE), referenced)
  if (length(fun_names) > 0) {
    keep_fun <- vapply(fun_names, function(v) {
      exists(v, envir = envir, inherits = FALSE) &&
        is.function(get(v, envir = envir))
    }, logical(1))
    user_funs <- fun_names[keep_fun]
  } else {
    user_funs <- character(0)
  }

  all_exports <- unique(c(export, auto_export, user_funs))

  if (nCores > 1) {
    # Create cluster:
    cl <- makePSOCKcluster(nCores)
    on.exit(stopCluster(cl), add = TRUE)

    # Reproducible parallel RNG streams:
    if (!is.null(seed)) {
      clusterSetRNGStream(cl, iseed = seed)
    }

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

    # Run with progress bar (the replication index 'i' is visible in expr):
    results <- pblapply(seq_len(reps), function(i) {
      eval(expr_sub)
    }, cl = cl)
  } else {
    # Reproducible sequential RNG:
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Sequential with progress bar. As in the parallel branch, the replication
    # index 'i' is visible in expr; all other symbols resolve in 'envir':
    results <- pblapply(seq_len(reps), function(i) {
      eval(expr_sub, envir = list2env(list(i = i), parent = envir))
    })
  }

  return(results)
}
