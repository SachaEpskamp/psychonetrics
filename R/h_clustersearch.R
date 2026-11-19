# Agglomerative search for clusters of variables.
#
# The search starts with the saturated model, in which every variable is its
# own construct and all constructs are freely related. In every step, each
# pair of constructs is combined into one cluster: a composite (cvm family) or
# a latent variable (lvm family) of all variables of the two constructs. The
# combination with the lowest value of the criterion becomes the new model, and
# the search continues with one construct less.
#
# All models are built with the regular model functions (cvm() or lvm()) and
# estimated with runmodel(). Only the models that are kept (one per step) are
# estimated in full; candidates are estimated without the baseline and
# saturated models, standard errors and modification indices, as only their
# likelihood and number of parameters are needed.

clustersearch <- function(
    data,
    vars,
    cluster = c("composite", "latent"),
    type = c("cov", "ggm"),
    criterion = "bic",
    min_clusters,
    improper = FALSE,
    nCores = 1,
    verbose = TRUE,
    ...
){
  experimentalWarning("clustersearch()")
  cluster <- match.arg(cluster)
  type <- match.arg(type)

  # The composite search uses the (experimental) cvm family. Show that note
  # here, as messages of the model functions are suppressed below:
  if (cluster == "composite"){
    experimentalWarning("cvm() / cnm() composite variable model family")
  }

  # Criterion (lower is better). Only criteria that are a function of the
  # likelihood and the number of parameters can be used:
  if (!is.character(criterion) || length(criterion) != 1){
    stop("'criterion' must be a single string.")
  }
  if (criterion == "aic") criterion <- "aic.ll"
  allowed_criteria <- c("bic", "aic.ll", "aic.ll2", "bic2", "ebic.25", "ebic.5", "ebic.75", "ebic1")
  if (!criterion %in% allowed_criteria){
    stop(paste0("'criterion' must be one of: ", paste0("'", allowed_criteria, "'", collapse = ", "), "."))
  }

  if (!is.logical(improper) || length(improper) != 1 || is.na(improper)){
    stop("'improper' must be TRUE or FALSE.")
  }
  if (!is.numeric(nCores) || length(nCores) != 1 || is.na(nCores) || nCores < 1){
    stop("'nCores' must be a positive number.")
  }

  # Arguments that are passed on to the model function:
  dots <- list(...)
  if (length(dots) > 0 && (is.null(names(dots)) || any(names(dots) == ""))){
    stop("All arguments in '...' must be named.")
  }
  dataArgNames <- c("groupvar", "groups", "covs", "means", "nobs", "missing", "covtype",
                    "standardize", "storedata", "estimator")
  modelArgNames <- c("optimizer", "identification", "equal")
  if (cluster == "composite") modelArgNames <- c(modelArgNames, "bad_marker_crit", "sumscores")
  notAllowed <- setdiff(names(dots), c(dataArgNames, modelArgNames))
  if (length(notAllowed) > 0){
    stop(paste0("The following argument(s) cannot be used in clustersearch(): ",
                paste0("'", notAllowed, "'", collapse = ", "),
                ". The search sets up the model itself; only arguments on the data, the estimator, ",
                "the optimizer and the identification are passed on to ",
                if (cluster == "composite") "cvm()." else "lvm()."))
  }
  if (!is.null(dots$estimator) && !dots$estimator %in% c("default", "ML", "FIML")){
    stop("clustersearch() only supports the estimators 'ML' and 'FIML'.")
  }
  dataargs <- dots[names(dots) %in% dataArgNames]
  modelargs <- dots[names(dots) %in% modelArgNames]

  # In the latent search the latent variables are by default identified by
  # fixing their variances to one, such that no cluster depends on the loading
  # of its first variable:
  if (cluster == "latent" && is.null(modelargs$identification)){
    modelargs$identification <- "variance"
  }

  # Variables:
  if (!missing(data)){
    if (is.matrix(data)) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("'data' must be a data frame or matrix.")
    dataargs$data <- data
    if (missing(vars)){
      if (is.null(colnames(data))) stop("'data' must have column names.")
      vars <- colnames(data)
      for (arg in c("groupvar", "groups")){
        if (is.character(dots[[arg]]) && length(dots[[arg]]) == 1) vars <- setdiff(vars, dots[[arg]])
      }
    }
  } else {
    if (is.null(dots$covs)) stop("Either 'data' or 'covs' (with 'nobs') must be supplied.")
    covs <- dots$covs
    covnames <- if (is.list(covs)) colnames(covs[[1]]) else colnames(covs)
    if (missing(vars)){
      if (is.null(covnames)) stop("'covs' must have column names, or 'vars' must be supplied.")
      vars <- covnames
    } else if (!is.null(covnames) && !identical(as.character(vars), covnames)){
      # Select the variables from the covariances and means:
      if (!all(vars %in% covnames)) stop("'vars' must be names of variables in 'covs'.")
      idx <- match(vars, covnames)
      dataargs$covs <- if (is.list(covs)){
        lapply(covs, function(x) x[idx, idx, drop = FALSE])
      } else if (length(dim(covs)) == 3){
        covs[idx, idx, , drop = FALSE]
      } else {
        covs[idx, idx, drop = FALSE]
      }
      means <- dots$means
      if (!is.null(means)){
        dataargs$means <- if (is.list(means)){
          lapply(means, function(x) x[idx])
        } else if (is.matrix(means)){
          means[idx, , drop = FALSE]
        } else {
          means[idx]
        }
      }
    }
  }
  if (!is.character(vars) || anyDuplicated(vars) > 0){
    stop("'vars' must be a character vector with unique variable names.")
  }
  dataargs$vars <- vars
  nVar <- length(vars)

  # Number of clusters at which the search stops:
  if (missing(min_clusters)){
    min_clusters <- if (cluster == "composite") 2 else 1
  }
  if (!is.numeric(min_clusters) || length(min_clusters) != 1 || is.na(min_clusters) ||
      min_clusters < 1 || min_clusters != round(min_clusters)){
    stop("'min_clusters' must be a positive whole number.")
  }
  if (cluster == "composite" && min_clusters < 2 && !isTRUE(modelargs$sumscores)){
    stop("'min_clusters' must be at least 2 in a search for composites: a single composite of all variables is the saturated model again, and its weights are not identified.")
  }
  if (nVar <= min_clusters){
    stop(paste0("The search needs more variables (", nVar, ") than 'min_clusters' (", min_clusters, ")."))
  }

  # Names of the clusters: the cluster formed in step s is named C<s>. Make
  # sure that these names do not clash with the names of the variables:
  prefix <- "C"
  while (any(paste0(prefix, seq_len(nVar)) %in% vars)) prefix <- paste0(prefix, "_")

  setup <- list(vars = vars, cluster = cluster, criterion = criterion,
                dataargs = dataargs, modelargs = modelargs,
                sampleStats = NULL, estimator = NULL)

  # Starting model: every variable is its own construct (saturated):
  partition <- stats::setNames(as.list(vars), vars)
  if (verbose) message("Estimating the starting model (", nVar, " variables)...")
  start <- clustersearch_fit(partition, setup, lean = FALSE)
  if (!is.na(start$error)){
    stop(paste0("The starting model could not be estimated: ", start$error))
  }
  if (!start$model@estimator %in% c("ML", "FIML")){
    stop("clustersearch() only supports the estimators 'ML' and 'FIML'.")
  }
  if (!is.finite(start$criterion)){
    stop(paste0("The criterion '", criterion, "' could not be computed for the starting model."))
  }

  # All further models use the sample statistics of the starting model:
  setup$sampleStats <- start$model@sample
  setup$estimator <- start$model@estimator
  setup$dataargs <- NULL

  # The starting model is the saturated model. Its log likelihood and number of
  # parameters give the chi-square and degrees of freedom of all other models:
  satLL <- start$logl
  satNpar <- start$npar

  # Multi-core:
  nCores <- round(nCores)
  cl <- NULL
  if (nCores > 1){
    cl <- parallel::makePSOCKcluster(nCores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library("psychonetrics"))
  }

  # Results:
  models <- list(clustersearch_settype(start$model, setup, type))
  partitions <- list(partition)
  warns <- list(start$warnings)
  path <- data.frame(step = 0, nclusters = nVar, cluster = NA_character_,
                     cluster1 = NA_character_, cluster2 = NA_character_, size = NA_integer_,
                     chisq = 0, df = 0,
                     npar = start$npar, criterion = start$criterion, stringsAsFactors = FALSE)
  candidates <- list()
  step <- 0

  while (length(partition) > min_clusters){
    step <- step + 1
    nClus <- length(partition)
    pairs <- utils::combn(nClus, 2)
    newname <- paste0(prefix, step)

    # All candidate partitions:
    candPartitions <- lapply(seq_len(ncol(pairs)), function(r){
      i <- pairs[1, r]
      j <- pairs[2, r]
      new <- list(vars[vars %in% c(partition[[i]], partition[[j]])])
      names(new) <- newname
      c(partition[-c(i, j)], new)
    })

    # Estimate the candidates:
    if (!is.null(cl)){
      fits <- parallel::parLapply(cl, candPartitions, clustersearch_fit, setup = setup, lean = TRUE)
    } else {
      fits <- lapply(candPartitions, clustersearch_fit, setup = setup, lean = TRUE)
    }

    cand <- data.frame(
      step = step,
      cluster1 = names(partition)[pairs[1, ]],
      cluster2 = names(partition)[pairs[2, ]],
      size = sapply(candPartitions, function(x) length(x[[newname]])),
      npar = sapply(fits, "[[", "npar"),
      logl = sapply(fits, "[[", "logl"),
      criterion = sapply(fits, "[[", "criterion"),
      converged = sapply(fits, "[[", "converged"),
      improper = sapply(fits, "[[", "improper"),
      error = sapply(fits, "[[", "error"),
      stringsAsFactors = FALSE
    )
    cand$df <- satNpar - cand$npar
    cand$chisq <- -2 * (cand$logl - satLL)
    cand$eligible <- is.na(cand$error) & is.finite(cand$criterion) & cand$converged %in% TRUE &
      (improper | cand$improper %in% FALSE)
    cand$selected <- FALSE

    # Select the best candidate, and estimate that model in full. If the full
    # model cannot be estimated, move on to the next best candidate:
    selected <- NULL
    full <- NULL
    for (r in order(cand$criterion)){
      if (!cand$eligible[r]) next
      full <- clustersearch_fit(candPartitions[[r]], setup, lean = FALSE)
      if (is.na(full$error) && is.finite(full$criterion) && isTRUE(full$converged) &&
          (improper || identical(full$improper, FALSE))){
        selected <- r
        break
      }
      # Record why the candidate is not selected after all:
      cand$eligible[r] <- FALSE
      cand$error[r] <- full$error
      cand$converged[r] <- isTRUE(full$converged)
      if (!is.na(full$improper)) cand$improper[r] <- full$improper
    }
    candidates[[step]] <- cand

    if (is.null(selected)){
      warning(paste0("None of the ", nrow(cand), " candidate models in step ", step,
                     " could be estimated properly (see the 'candidates' element of the output). ",
                     "The search was stopped at ", nClus, " clusters."), call. = FALSE)
      step <- step - 1
      break
    }
    candidates[[step]]$selected[selected] <- TRUE

    partition <- candPartitions[[selected]]
    partitions[[step + 1]] <- partition
    models[[step + 1]] <- clustersearch_settype(full$model, setup, type)
    warns[[step + 1]] <- full$warnings
    path <- rbind(path, data.frame(step = step, nclusters = length(partition), cluster = newname,
                                   cluster1 = cand$cluster1[selected], cluster2 = cand$cluster2[selected],
                                   size = cand$size[selected],
                                   chisq = -2 * (full$logl - satLL), df = satNpar - full$npar,
                                   npar = full$npar, criterion = full$criterion, stringsAsFactors = FALSE))

    if (verbose){
      message("Step ", step, " (", length(partition), if (length(partition) == 1) " cluster): " else " clusters): ",
              cand$cluster1[selected], " + ",
              cand$cluster2[selected], " -> ", newname, " (", criterion, " = ",
              format(round(full$criterion, 2), nsmall = 2), "; ", sum(!cand$eligible),
              " of ", nrow(cand), " candidates skipped)")
    }
  }

  candidates <- do.call(rbind, candidates)
  if (!is.null(candidates)){
    candidates <- candidates[, c("step", "cluster1", "cluster2", "size", "chisq", "df", "npar",
                                 "criterion", "converged", "improper", "error", "eligible", "selected")]
    names(candidates)[names(candidates) == "criterion"] <- criterion
    rownames(candidates) <- NULL
  }
  names(path)[names(path) == "criterion"] <- criterion
  rownames(path) <- NULL
  names(models) <- names(partitions) <- names(warns) <- as.character(path$nclusters)

  best <- which.min(path[[criterion]])
  if (length(warns[[best]]) > 0){
    warning(paste0("The selected model gave the following warning(s) when it was estimated: ",
                   paste(warns[[best]], collapse = " | ")), call. = FALSE)
  }

  out <- list(
    model = models[[best]],
    models = models,
    path = path,
    candidates = candidates,
    clusters = partitions[[best]],
    partitions = partitions,
    warnings = warns,
    best = path$step[best],
    criterion = criterion,
    cluster = cluster,
    type = type,
    vars = vars
  )
  class(out) <- "clustersearch"
  out
}

# Build the model of a partition. 'partition' is a named list with the
# variables of every construct.
clustersearch_build <- function(partition, setup, lean = TRUE){
  vars <- setup$vars
  nVar <- length(vars)
  nClus <- length(partition)
  mat <- matrix(0, nVar, nClus)
  for (k in seq_len(nClus)) mat[match(partition[[k]], vars), k] <- 1

  args <- c(list(latents = names(partition), baseline_saturated = !lean, verbose = FALSE),
            setup$modelargs)
  if (is.null(setup$sampleStats)){
    args <- c(args, setup$dataargs)
  } else {
    args <- c(args, list(sampleStats = setup$sampleStats, estimator = setup$estimator))
  }

  if (setup$cluster == "composite"){
    # Composites: a variable that is not clustered is a composite with a
    # single indicator.
    do.call(cvm, c(list(weights = mat, latent = "cov", composite = "cov"), args))
  } else {
    # Latent variables: a variable that is not clustered is a latent variable
    # with a single indicator and no residual variance.
    clustered <- vars %in% unlist(partition[lengths(partition) > 1])
    eps <- diag(as.numeric(clustered), nVar, nVar)
    mod <- do.call(lvm, c(list(lambda = mat, latent = "cov", residual = "cov", sigma_epsilon = eps,
                               simplelambdastart = TRUE), args))
    clustersearch_latent_start(mod, partition, vars)
  }
}

# Standardized loadings of a single latent variable on the variables 'idx',
# from the correlation matrix R. Used for starting values only.
clustersearch_onefactor <- function(R, idx){
  n <- length(idx)
  Rk <- R[idx, idx, drop = FALSE]
  if (n == 2){
    # Two variables only determine the product of the loadings. Their relative
    # size follows from the correlations with all other variables:
    r <- Rk[1, 2]
    rest <- setdiff(seq_len(nrow(R)), idx)
    ratio <- 1
    signs <- c(1, sign(r))
    if (length(rest) > 0){
      u <- tryCatch(svd(R[idx, rest, drop = FALSE], nu = 1, nv = 0)$u[, 1], error = function(e) c(NA, NA))
      if (all(is.finite(u)) && u[1] * u[2] * r > 0){
        ratio <- min(max(abs(u[1] / u[2]), 1/3), 3)
      }
    }
    l <- signs * c(sqrt(abs(r) * ratio), sqrt(abs(r) / ratio))
  } else {
    # Principal axis factoring with one factor:
    h <- tryCatch(1 - 1 / diag(solve(Rk)), error = function(e) rep(NA, n))
    if (any(!is.finite(h))) h <- apply(abs(Rk - diag(n)), 1, max)
    h <- pmin(pmax(h, 0.05), 0.95)
    l <- sqrt(h)
    for (it in 1:50){
      diag(Rk) <- h
      e <- eigen(Rk, symmetric = TRUE)
      l <- e$vectors[, 1] * sqrt(max(e$values[1], 0))
      hNew <- pmin(pmax(l^2, 0.05), 0.95)
      if (max(abs(hNew - h)) < 1e-4) break
      h <- hNew
    }
    if (sum(l) < 0) l <- -l
  }
  s <- sign(l)
  s[s == 0] <- 1
  s * pmin(pmax(abs(l), 0.1), 0.95)
}

# Starting values for the latent variable model of a partition.
#
# The starting values of lvm() come from an exploratory factor analysis with as
# many factors as there are latent variables. In this search most latent
# variables have a single indicator, for which such an analysis is slow and
# gives poor values. The model is therefore built with simple starting values,
# which are replaced here: one-factor loadings per cluster, the correlations
# between the constructs that follow from these loadings, and the remaining
# variances as residual variances.
clustersearch_latent_start <- function(x, partition, vars){
  pars <- x@parameters
  nVar <- length(vars)
  nClus <- length(partition)
  member <- lapply(partition, match, table = vars)
  size <- lengths(member)

  for (g in seq_len(nrow(x@sample@groups))){
    S <- as.matrix(x@sample@covs[[g]])
    if (any(!is.finite(S)) || any(diag(S) <= 0)) next
    S <- as.matrix(spectralshift(S))
    sds <- sqrt(diag(S))
    R <- stats::cov2cor(S)

    # Standardized loadings:
    l <- rep(1, nVar)
    for (k in which(size > 1)) l[member[[k]]] <- clustersearch_onefactor(R, member[[k]])

    # Correlations between the constructs:
    Phi <- diag(nClus)
    if (nClus > 1){
      for (k in seq_len(nClus - 1)){
        for (j in (k + 1):nClus){
          a <- l[member[[k]]]
          b <- l[member[[j]]]
          Phi[k, j] <- Phi[j, k] <- sum(outer(a, b) * R[member[[k]], member[[j]], drop = FALSE]) /
            (sum(a^2) * sum(b^2))
        }
      }
      off <- row(Phi) != col(Phi)
      Phi[off] <- pmin(pmax(Phi[off], -0.95), 0.95)
      it <- 0
      while (min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) < 0.05 && it < 100){
        Phi[off] <- 0.9 * Phi[off]
        it <- it + 1
      }
    }

    # Loadings, variances and residual variances:
    Lambda <- matrix(0, nVar, nClus)
    for (k in seq_len(nClus)) Lambda[member[[k]], k] <- l[member[[k]]] * sds[member[[k]]]
    theta <- diag(S) * (1 - l^2)

    # Scale to the identification of the model: a latent variable with a fixed
    # (nonzero) loading is scaled by that loading.
    inGroup <- pars$group_id == g
    scale <- rep(1, nClus)
    for (k in seq_len(nClus)){
      marker <- which(inGroup & pars$matrix == "lambda" & pars$col == k & pars$fixed & pars$est != 0)
      if (length(marker) > 0){
        i <- pars$row[marker[1]]
        sc <- Lambda[i, k] / pars$est[marker[1]]
        if (abs(sc) < 0.05 * sds[i]) sc <- 0.05 * sds[i] * ifelse(sc < 0, -1, 1)
        scale[k] <- sc
      }
    }
    Lambda <- Lambda %*% diag(1 / scale, nClus, nClus)
    Psi <- diag(scale, nClus, nClus) %*% Phi %*% diag(scale, nClus, nClus)

    rows <- which(inGroup & pars$matrix == "lambda" & !pars$fixed)
    pars$est[rows] <- Lambda[cbind(pars$row[rows], pars$col[rows])]
    rows <- which(inGroup & pars$matrix == "sigma_zeta" & !pars$fixed)
    pars$est[rows] <- Psi[cbind(pars$row[rows], pars$col[rows])]
    rows <- which(inGroup & pars$matrix == "sigma_epsilon" & !pars$fixed & pars$row == pars$col)
    pars$est[rows] <- theta[pars$row[rows]]
  }
  x@parameters <- pars

  # Make the estimates of equality constrained parameters consistent:
  updateModel(parVector(x), x, updateMatrices = FALSE)
}

# Did the optimizer converge? Uses the convergence code of the optimizer and,
# with gradient = TRUE, the same gradient criterion as runmodel(). Evaluating
# the gradient takes about as long as an iteration of the optimizer, so it is
# only checked for the models that are kept.
clustersearch_converged <- function(x, gradient = TRUE){
  conv <- x@optim$convergence
  if (!is.null(conv) && length(conv) == 1 && !is.na(conv) && is.numeric(conv) && conv != 0){
    return(FALSE)
  }
  if (!gradient) return(TRUE)
  grad <- tryCatch({
    if (x@cpp) psychonetrics_gradient_cpp(parVector(x), x) else psychonetrics_gradient(parVector(x), x)
  }, error = function(e) NA)
  if (x@cpp) clearPrepCache_cpp()
  if (any(!is.finite(grad))) return(FALSE)
  mean(abs(grad)) <= 1
}

# Is the solution improper? Only relevant for latent variables: a parameter on
# its bound (such as a residual variance of zero) or a variance-covariance
# matrix of the constructs that is not positive definite.
clustersearch_improper <- function(x, cluster){
  if (cluster != "latent") return(FALSE)
  pars <- x@parameters
  free <- !pars$fixed
  tol <- sqrt(.Machine$double.eps)
  if (any(pars$est[free] <= pars$minimum[free] + tol) || any(pars$est[free] >= pars$maximum[free] - tol)){
    return(TRUE)
  }
  for (g in seq_along(x@modelmatrices)){
    Psi <- as.matrix(x@modelmatrices[[g]]$sigma_zeta)
    if (any(!is.finite(Psi)) || any(diag(Psi) <= 0)) return(TRUE)
    ev <- eigen(stats::cov2cor(Psi), symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= tol) return(TRUE)
    Theta <- diag(as.matrix(x@modelmatrices[[g]]$sigma_epsilon))
    if (any(!is.finite(Theta)) || any(Theta < 0)) return(TRUE)
  }
  FALSE
}

# Build and estimate the model of a partition. With lean = TRUE only the
# likelihood is obtained (no baseline and saturated model, standard errors or
# modification indices), and the model itself is not returned.
clustersearch_fit <- function(partition, setup, lean = TRUE){
  warns <- character(0)
  res <- withCallingHandlers(
    tryCatch({
      mod <- suppressMessages(clustersearch_build(partition, setup, lean = lean))
      if (lean){
        suppressMessages(runmodel(mod, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE,
                                  log = FALSE, warn_gradient = FALSE, warn_bounds = FALSE))
      } else {
        suppressMessages(runmodel(mod, warn_gradient = FALSE, warn_bounds = FALSE))
      }
    }, error = function(e) e),
    warning = function(w){
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  out <- list(model = NULL, error = NA_character_, warnings = unique(warns), logl = NA_real_,
              npar = NA_real_, criterion = NA_real_, converged = NA, improper = NA)
  if (inherits(res, "error")){
    out$error <- conditionMessage(res)
    return(out)
  }
  fm <- res@fitmeasures
  if (!is.null(fm$logl)) out$logl <- fm$logl
  if (!is.null(fm$npar)) out$npar <- fm$npar
  if (!is.null(fm[[setup$criterion]])) out$criterion <- fm[[setup$criterion]]
  out$converged <- clustersearch_converged(res, gradient = !lean)
  out$improper <- clustersearch_improper(res, setup$cluster)
  if (!lean) out$model <- res
  out
}

# Transform an estimated model to the requested type. All models are estimated
# with variance-covariance matrices (type = "cov"); for type = "ggm" the
# equivalent network model is returned (cnm / lnm).
clustersearch_settype <- function(x, setup, type){
  if (type == "cov") return(x)
  res <- tryCatch(suppressWarnings(suppressMessages({
    if (setup$cluster == "composite"){
      mod <- transmod(x, latent = "ggm", composite = "ggm", verbose = FALSE)
    } else {
      mod <- transmod(x, latent = "ggm", verbose = FALSE)
    }
    runmodel(mod)
  })), error = function(e) e)
  if (inherits(res, "error")){
    warning(paste0("A model could not be transformed to type = 'ggm' and is returned with type = 'cov': ",
                   conditionMessage(res)), call. = FALSE)
    return(x)
  }
  res
}

# Print method:
print.clustersearch <- function(x, digits = 2, ...){
  crit <- x$criterion
  path <- x$path
  nVar <- length(x$vars)

  cat("Cluster search (psychonetrics)\n")
  cat("  Clusters:", if (x$cluster == "composite") "composites" else "latent variables",
      if (x$cluster == "composite"){
        if (x$type == "ggm") "(cnm)" else "(cvm)"
      } else {
        if (x$type == "ggm") "(lnm)" else "(lvm)"
      }, "\n")
  cat("  Variables:", nVar, "\n")
  cat("  Criterion:", crit, "\n\n")

  merged <- ifelse(is.na(path$cluster), "(all variables separate)",
                   paste0(path$cluster1, " + ", path$cluster2, " -> ", path$cluster))
  tab <- data.frame(
    step = path$step,
    clusters = path$nclusters,
    merged = format(c("merged", merged))[-1],
    chisq = formatC(path$chisq, digits = digits, format = "f"),
    df = path$df,
    crit = formatC(path[[crit]], digits = digits, format = "f"),
    best = ifelse(path$step == x$best, "<-", "  "),
    stringsAsFactors = FALSE
  )
  names(tab)[names(tab) == "merged"] <- format(c("merged", merged))[1]
  names(tab)[names(tab) == "crit"] <- crit
  names(tab)[names(tab) == "best"] <- "  "
  print(tab, row.names = FALSE)

  cat("\nSelected model: step ", x$best, " (", length(x$clusters),
      if (length(x$clusters) == 1) " cluster, " else " clusters, ", crit, " = ",
      formatC(min(path[[crit]]), digits = digits, format = "f"), ")\n", sep = "")
  multi <- lengths(x$clusters) > 1
  if (any(multi)){
    for (k in which(multi)){
      cat("  ", names(x$clusters)[k], ": ", paste(x$clusters[[k]], collapse = ", "), "\n", sep = "")
    }
  }
  if (any(!multi)){
    cat("  Not clustered: ", paste(unlist(x$clusters[!multi]), collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$candidates)){
    nSkip <- sum(!x$candidates$eligible)
    if (nSkip > 0){
      cat("\n", nSkip, " of ", nrow(x$candidates),
          " candidate models were skipped (not estimated properly", if (x$cluster == "latent") " or improper",
          "); see the 'candidates' element.\n", sep = "")
    }
  }
  cat("\nUse $model for the selected model, $models for the model of every step, and plot() for the dendrogram.\n")
  invisible(x)
}

# The tree of a search: one node per variable and per cluster, with the step
# in which the node is formed, its children and its vertical position.
clustersearch_tree <- function(x){
  vars <- x$vars
  path <- x$path[x$path$step > 0, , drop = FALSE]
  nodes <- data.frame(name = vars, step = 0, child1 = NA_integer_, child2 = NA_integer_,
                      stringsAsFactors = FALSE)
  for (s in seq_len(nrow(path))){
    nodes <- rbind(nodes, data.frame(name = path$cluster[s], step = path$step[s],
                                     child1 = match(path$cluster1[s], nodes$name),
                                     child2 = match(path$cluster2[s], nodes$name),
                                     stringsAsFactors = FALSE))
  }
  # Nodes that are not merged into another node:
  top <- setdiff(seq_len(nrow(nodes)), c(nodes$child1, nodes$child2))
  # (within a cluster, the part that contains the earliest variable comes first)
  leaves <- function(i){
    if (is.na(nodes$child1[i])) return(i)
    a <- leaves(nodes$child1[i])
    b <- leaves(nodes$child2[i])
    if (min(a) < min(b)) c(a, b) else c(b, a)
  }
  # Order the top nodes by their first variable:
  firstLeaf <- sapply(top, function(i) min(leaves(i)))
  top <- top[order(firstLeaf)]
  ord <- unlist(lapply(top, leaves))
  nodes$y <- NA_real_
  nodes$y[ord] <- rev(seq_along(ord))
  for (i in which(!is.na(nodes$child1))){
    nodes$y[i] <- mean(nodes$y[c(nodes$child1[i], nodes$child2[i])])
  }
  nodes$top <- seq_len(nrow(nodes)) %in% top
  nodes
}

# Plot method: dendrogram with the criterion on the same horizontal axis.
plot.clustersearch <- function(x, labels, cex = 1, nodelabels = TRUE, col = "black",
                               grid = "grey90", ...){
  vars <- x$vars
  nVar <- length(vars)
  path <- x$path
  crit <- x$criterion
  nSteps <- max(path$step)
  if (nSteps < 1) stop("The search did not make any step: nothing to plot.")
  if (missing(labels)) labels <- vars
  if (length(labels) != nVar) stop("'labels' must have one element per variable.")
  nodes <- clustersearch_tree(x)

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(op)
  })
  graphics::layout(matrix(1:2, 2, 1), heights = c(max(1.25, nVar / 6), 1))

  # Same left and right margins in both panels, wide enough for the labels:
  left <- max(graphics::strwidth(labels, units = "inches", cex = 0.9 * cex)) + 0.35
  left <- max(left, 1)
  right <- 0.25

  ### Dendrogram ###
  graphics::par(mai = c(0.1, left, 0.3, right))
  graphics::plot(NA, xlim = c(0, nSteps), ylim = c(0.5, nVar + 0.5), axes = FALSE, xlab = "", ylab = "",
                 xaxs = "r", yaxs = "i")
  graphics::abline(v = 0:nSteps, col = grid)
  graphics::abline(v = x$best, col = col, lty = 2)
  for (i in seq_len(nrow(nodes))){
    kids <- c(nodes$child1[i], nodes$child2[i])
    if (!is.na(kids[1])){
      graphics::segments(nodes$step[kids], nodes$y[kids], nodes$step[i], nodes$y[kids], col = col, lwd = 1.5)
      graphics::segments(nodes$step[i], nodes$y[kids[1]], nodes$step[i], nodes$y[kids[2]], col = col, lwd = 1.5)
    }
    # Nodes that are not merged run to the end of the search:
    if (nodes$top[i] && nodes$step[i] < nSteps){
      graphics::segments(nodes$step[i], nodes$y[i], nSteps, nodes$y[i], col = col, lwd = 1.5)
    }
  }
  isCluster <- !is.na(nodes$child1)
  graphics::points(nodes$step[isCluster], nodes$y[isCluster], pch = 21, bg = col, col = "white", cex = 1.4 * cex)
  if (nodelabels){
    graphics::text(nodes$step[isCluster], nodes$y[isCluster], labels = nodes$name[isCluster],
                   pos = 2, offset = 0.45, cex = 0.75 * cex, col = col, xpd = NA)
  }
  leaf <- which(nodes$step == 0)
  graphics::axis(2, at = nodes$y[leaf], labels = labels, las = 1, tick = FALSE, line = -0.6, cex.axis = 0.9 * cex)

  ### Criterion ###
  graphics::par(mai = c(0.9, left, 0.1, right))
  yr <- range(path[[crit]], finite = TRUE)
  if (diff(yr) == 0) yr <- yr + c(-1, 1)
  ticks <- pretty(yr)
  graphics::plot(NA, xlim = c(0, nSteps), ylim = range(c(yr, ticks)), axes = FALSE,
                 xlab = "Number of clusters", ylab = "", xaxs = "r", cex.lab = cex)
  graphics::abline(v = 0:nSteps, col = grid)
  graphics::abline(h = ticks, col = grid)
  graphics::abline(v = x$best, col = col, lty = 2)
  graphics::lines(path$step, path[[crit]], col = col, lwd = 2)
  graphics::points(path$step, path[[crit]], pch = 21, bg = col, col = "white", cex = 1.2 * cex)
  graphics::points(x$best, min(path[[crit]]), pch = 21, bg = "white", col = col, cex = 1.9 * cex, lwd = 2)
  graphics::axis(1, at = path$step, labels = path$nclusters, cex.axis = 0.9 * cex)
  graphics::axis(2, at = ticks, las = 1, cex.axis = 0.9 * cex)
  graphics::mtext(toupper(crit), side = 2, line = 3.6, cex = cex)

  invisible(nodes)
}
