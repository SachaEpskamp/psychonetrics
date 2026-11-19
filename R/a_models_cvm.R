# Composite variable model (cvm): structural equation models with latent
# variables (common factors) and composites, after Schamberger, Schuberth,
# Henseler & Rosseel (factor- and composite-based SEM).

# Pattern of a typed covariance block of the composite indicators. Composite
# indicators of different composites cannot be connected in these matrices
# (their covariances are implied by the relations between the composites).
cvm_blockpattern <- function(x, block, nGroup, name, kind = c("symmetric", "offdiagonal", "diagonal")){
  kind <- match.arg(kind)
  n <- length(block)
  same <- outer(block, block, "==")

  if (is.character(x)){
    x <- x[1]
    if (x == "empty") x <- "zero"
    if (!x %in% c("full", "diag", "zero")){
      stop(paste0("'", name, "' must be 'full', 'diag', 'zero' or a matrix."))
    }
    mat <- matrix(0, n, n)
    if (kind == "diagonal"){
      if (x != "zero") diag(mat) <- 1
    } else if (kind == "offdiagonal"){
      if (x == "full") mat[same] <- 1
      diag(mat) <- 0
    } else {
      if (x == "full") mat[same] <- 1
      if (x == "diag") diag(mat) <- 1
    }
    return(mat)
  }

  # Matrix or array input:
  if (!is.numeric(x) && !is.logical(x)) stop(paste0("'", name, "' must be 'full', 'diag', 'zero' or a matrix."))
  x <- 1 * x
  if (length(dim(x)) == 2){
    x <- array(x, c(dim(x), nGroup))
  }
  if (dim(x)[1] != n || dim(x)[2] != n){
    stop(paste0("'", name, "' must be a square matrix with one row and column per composite indicator (", n, ")."))
  }
  warned <- FALSE
  for (g in seq_len(dim(x)[3])){
    cur <- x[,,g]
    if (any(cur[!same] != 0)){
      if (!warned){
        warning(paste0("Elements of '", name, "' connecting indicators of different composites are fixed to zero: ",
                       "these (partial) covariances are implied by the relations between the composites."))
        warned <- TRUE
      }
      cur[!same] <- 0
      x[,,g] <- cur
    }
  }
  x
}

# OLS based starting values for beta, given the (total) covariances of the constructs:
cvm_beta_start <- function(betaFixed, expLatSigma){
  nGroup <- dim(betaFixed)[3]
  nLatent <- dim(betaFixed)[1]
  beta_start_list <- vector("list", nGroup)

  for (g in seq_len(nGroup)){
    S <- as.matrix(expLatSigma[[g]])
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || min(ev) / max(ev) < 1e-8) return(NULL)

    betaMat <- betaFixed[,,g]
    beta_g <- matrix(0, nLatent, nLatent)
    for (i in seq_len(nLatent)){
      predictors <- which(betaMat[i,] != 0)
      if (length(predictors) == 0) next
      S_pp <- S[predictors, predictors, drop = FALSE]
      s_pi <- S[predictors, i, drop = FALSE]
      if (length(predictors) > 1){
        ev_pp <- eigen(S_pp, symmetric = TRUE, only.values = TRUE)$values
        if (any(!is.finite(ev_pp)) || max(ev_pp) / min(ev_pp) > 1e6) return(NULL)
      }
      ols <- tryCatch(as.numeric(solve(S_pp, s_pi)), error = function(e) NULL)
      if (is.null(ols) || any(!is.finite(ols))) return(NULL)
      beta_g[i, predictors] <- ols
    }
    if (any(abs(beta_g[betaMat != 0]) > 5)) return(NULL)
    IminB_inv <- tryCatch(solve(diag(nLatent) - beta_g), error = function(e) NULL)
    if (is.null(IminB_inv) || max(abs(IminB_inv)) > 100) return(NULL)
    beta_start_list[[g]] <- beta_g
  }
  beta_start_list
}

# Starting values for the composite weights, per group: partial least squares
# with mode B outer estimation and the factorial inner scheme. Every composite
# is repeatedly regressed on the (proxies of the) constructs it is related to
# in the model. For two composites this is canonical correlation analysis; in
# general it gives consistent estimates of the weights of the composite model
# (Dijkstra, 2017), and therefore starting values close to the optimum.
#   S:         covariance matrix of the observed variables
#   block:     composite (column) of every variable, 0 if not a composite indicator
#   comp:      logical, which constructs are composites
#   facProxy:  variables x constructs matrix with proxy weights for the factors
#   adjacency: constructs x constructs logical matrix, which constructs are related
#   genvar:    logical, is the construct part of the model saturated? The mode B
#              weights are then improved with cvm_genvar_weights()
# Returns a variables x constructs matrix in which the column of every
# composite holds its weights, scaled to a composite with unit variance.
cvm_start_weights <- function(S, block, comp, facProxy, adjacency, maxit = 200, tol = 1e-8, genvar = FALSE){
  nLat <- length(comp)
  nVar <- nrow(S)
  sds <- sqrt(diag(S))
  compIdx <- which(comp)

  # Scale a vector of weights to a proxy with unit variance:
  stdize <- function(v){
    s <- drop(t(v) %*% S %*% v)
    if (is.finite(s) && s > 0) v / sqrt(s) else v
  }

  # Start with equal standardized weights:
  V <- facProxy
  for (k in compIdx){
    idx <- which(block == k)
    V[idx, k] <- 1 / sds[idx]
  }
  for (j in seq_len(nLat)) V[, j] <- stdize(V[, j])
  Vunit <- V

  # Composites with more than one indicator:
  multi <- compIdx[vapply(compIdx, function(k) sum(block == k) > 1, logical(1))]
  ok <- TRUE
  if (length(multi) > 0){
    for (it in seq_len(maxit)){
      Vold <- V
      for (k in multi){
        idx <- which(block == k)
        nb <- setdiff(which(adjacency[k, ]), k)
        if (length(nb) == 0) nb <- setdiff(seq_len(nLat), k)
        if (length(nb) == 0) next
        # Inner weights (factorial scheme) and the inner proxy:
        e <- drop(t(V[, k]) %*% S %*% V[, nb, drop = FALSE])
        z <- V[, nb, drop = FALSE] %*% e
        # Outer weights (mode B): regression of the inner proxy on the indicators:
        w <- tryCatch(as.vector(solve(S[idx, idx, drop = FALSE], S[idx, , drop = FALSE] %*% z)),
                      error = function(e) NULL)
        if (is.null(w) || any(!is.finite(w)) || all(abs(w) < .Machine$double.eps)){
          ok <- FALSE
          break
        }
        v <- numeric(nVar)
        v[idx] <- w
        V[, k] <- stdize(v)
      }
      if (!ok || any(!is.finite(V))){
        ok <- FALSE
        break
      }
      if (max(abs(abs(V) - abs(Vold))) < tol) break
    }
  }

  # Fall back to equal standardized weights:
  if (!ok) V <- Vunit

  # Saturated construct part: the maximum likelihood weights have a closed form
  # criterion, see cvm_genvar_weights():
  if (genvar && length(multi) > 0){
    V <- cvm_genvar_weights(S, block, V, multi)
  }
  V
}

# Weights that minimize the determinant of the correlation matrix of the
# composites (the GENVAR criterion of Kettenring, 1971, Biometrika), used as
# starting values when the construct part of the model is saturated: all
# constructs are composites that are all related to each other (covariances or
# a network without zeros, and no effects in beta).
#
# In that model, with saturated covariances within the composites, the
# chi-square is
#   N * (sum_k log|S_kk| + log|R_C| - log|S|),
# with S_kk the sample covariances of the indicators of composite k and R_C the
# correlation matrix of the composites. The maximum likelihood estimates of the
# weights therefore minimize |R_C| (and the covariances within and between the
# composites are S_kk and W'SW). The mode B weights are consistent, but do not
# minimize this criterion; starting from them, the optimizer needed many
# iterations in large models (114 iterations and about an hour for 83 variables
# in 9 composites, against 8 iterations from the weights below).
#
# The weights are improved one composite at a time. Given the other
# composites, |R_C| is minimal when the multiple correlation of the composite
# with the other composites is maximal: the first canonical variate of its
# indicators with the other composites. Every step lowers the criterion.
#   V:     variables x constructs matrix of starting weights (unit variance)
#   multi: the composites with more than one indicator
# Returns V with improved weights (unit variance), or V itself if the
# criterion could not be improved.
cvm_genvar_weights <- function(S, block, V, multi, maxit = 1000, tol = 1e-10){
  nLat <- ncol(V)
  if (nLat < 2 || length(multi) == 0) return(V)

  # Log determinant of the correlation matrix of the constructs:
  criterion <- function(V){
    R <- cov2cor(t(V) %*% S %*% V)
    d <- determinant(R, logarithm = TRUE)
    if (d$sign > 0) as.numeric(d$modulus) else NA_real_
  }

  Vstart <- V
  critStart <- tryCatch(criterion(V), error = function(e) NA_real_)
  if (!is.finite(critStart)) return(Vstart)

  res <- tryCatch({
    # Cholesky factors of the covariances within the composites:
    idxList <- lapply(seq_len(nLat), function(k) which(block == k))
    cholList <- lapply(seq_len(nLat), function(k){
      if (k %in% multi) chol(S[idxList[[k]], idxList[[k]], drop = FALSE]) else NULL
    })
    SV <- S %*% V
    critOld <- critStart

    for (it in seq_len(maxit)){
      for (k in multi){
        idx <- idxList[[k]]
        others <- setdiff(seq_len(nLat), k)
        # Covariances of the indicators with the other constructs, and the
        # covariances of the other constructs:
        Cko <- SV[idx, others, drop = FALSE]
        Voo <- crossprod(V[, others, drop = FALSE], SV[, others, drop = FALSE])
        A <- Cko %*% solve(Voo, t(Cko))
        # First canonical variate: largest eigenvector of S_kk^-1 A:
        R <- cholList[[k]]
        M <- backsolve(R, t(backsolve(R, A, transpose = TRUE)), transpose = TRUE)
        u <- eigen(0.5 * (M + t(M)), symmetric = TRUE)$vectors[, 1]
        w <- backsolve(R, u)
        # Keep the orientation of the composite:
        if (sum(w * SV[idx, k]) < 0) w <- -w
        V[idx, k] <- w
        SV[, k] <- S[, idx, drop = FALSE] %*% w
      }
      critNew <- criterion(V)
      if (!is.finite(critNew)) stop("criterion not finite")
      converged <- critOld - critNew < tol
      critOld <- critNew
      if (converged) break
    }
    list(V = V, crit = critOld)
  }, error = function(e) NULL)

  if (is.null(res) || any(!is.finite(res$V)) || !(res$crit <= critStart)) return(Vstart)
  res$V
}

# Check the indicators of every composite for linear dependencies, using the
# eigenvalues of their sample correlation matrix:
#  - (numerically) singular: the indicators are linearly dependent (for example
#    a total score next to its items, or more indicators than observations).
#    The weights of such a composite are not identified and the covariances
#    within the composite cannot be estimated, so an error is given;
#  - condition number above 'crit': the indicators are nearly linearly
#    dependent. The model can be estimated, but the weights are poorly
#    determined (as regression weights of collinear predictors), so a warning
#    is given.
#   covs:  list with the sample covariance matrix of every group
#   block: composite (column) of every variable, 0 if not a composite indicator
cvm_check_collinearity <- function(covs, block, compIdx, latents, sumscores = FALSE, crit = 1000){
  dependent <- character(0)
  cond <- numeric(0)
  for (g in seq_along(covs)){
    S <- as.matrix(covs[[g]])
    for (k in compIdx){
      idx <- which(block == k)
      if (length(idx) < 2) next
      Sk <- S[idx, idx, drop = FALSE]
      if (any(!is.finite(Sk)) || any(diag(Sk) <= 0)) next
      ev <- eigen(cov2cor(Sk), symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) <= sqrt(.Machine$double.eps) * max(ev)){
        dependent <- c(dependent, latents[k])
      } else if (max(ev) / min(ev) > crit){
        cond[latents[k]] <- max(cond[latents[k]], max(ev) / min(ev), na.rm = TRUE)
      }
    }
  }
  if (length(dependent) > 0){
    stop(paste0("The indicators of the following composite(s) are linearly dependent (their sample correlation matrix is singular): ",
                paste(unique(dependent), collapse = ", "),
                ". A composite variable model cannot be estimated for such a composite. Remove indicators that are a combination of other ",
                "indicators of the same composite (for example a total score next to its items), and make sure that a composite ",
                "does not have more indicators than there are observations."), call. = FALSE)
  }
  if (length(cond) > 0){
    warning(paste0("The indicators of the following composite(s) are nearly linearly dependent (condition number of their sample correlation matrix in parentheses): ",
                   paste0(names(cond), " (", format(signif(cond, 3), scientific = FALSE, trim = TRUE), ")", collapse = ", "), ". ",
                   if (!sumscores) "The weights of such a composite are poorly determined, and the " else "The ",
                   "optimizer may not report convergence."), call. = FALSE)
  }
  invisible(NULL)
}

# Starting values of the covariance structure of the composite indicators (the
# "_c" matrices) for the types "ggm", "prec" and "chol", formed per composite.
#
# The "_c" matrices are block diagonal, with one block per composite. The
# matrixsetup functions that are shared with the other families form starting
# values for all variables at once: with glasso and zeros between the
# composites for the network and precision types, and from the Cholesky factor
# of the full matrix for the Cholesky type. With highly correlated indicators
# within a composite glasso can then return a matrix that is not positive
# definite (the starting values were improper, after which estimation started
# from the emergency starting values of runmodel() and stopped far from the
# optimum in examples), and the blocks of the Cholesky factor of a matrix are
# not the Cholesky factors of its blocks. Here every block is treated on its
# own:
#  - no zeros within the block: the inverse (or Cholesky factor) of the sample
#    covariances of the indicators of the composite;
#  - zeros within the block (network and precision types): glasso on the block;
#  - if the result is not positive definite: independent indicators.
#   expcov: list with the sample covariances of the composite indicators per group
#   block:  composite of every composite indicator
cvm_block_start <- function(modMatrices, type, expcov, block, nGroup){
  if (!type %in% c("ggm", "prec", "chol")) return(modMatrices)

  # Precision matrix of a block, given which off-diagonal elements are in the model:
  blockPrecision <- function(S, free){
    p <- nrow(S)
    independent <- diag(1 / diag(S), p, p)
    if (p == 1) return(independent)
    zeroes <- which(!free & diag(p) != 1, arr.ind = TRUE)
    K <- tryCatch({
      if (nrow(zeroes) == 0){
        solve(S)
      } else if (nrow(zeroes) == p * (p - 1)){
        independent
      } else {
        glasso(S, rho = 1e-10, zero = zeroes)$wi
      }
    }, error = function(e) NULL)
    if (is.null(K) || any(!is.finite(K))) return(independent)
    K <- 0.5 * (K + t(K))
    K[!free & diag(p) != 1] <- 0
    ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= sqrt(.Machine$double.eps) * max(abs(ev))) return(independent)
    K
  }

  for (g in seq_len(nGroup)){
    S <- as.matrix(expcov[[g]])
    for (k in unique(block)){
      idx <- which(block == k)
      p <- length(idx)
      Sk <- S[idx, idx, drop = FALSE]
      if (any(!is.finite(Sk)) || any(diag(Sk) <= 0)) next

      if (type == "chol"){
        pattern <- matrix(modMatrices$lowertri_c[[1]][idx, idx, g] != 0, p, p)
        L <- tryCatch(t(chol(Sk)), error = function(e) diag(sqrt(diag(Sk)), p, p))
        modMatrices$lowertri_c$start[idx, idx, g] <- pattern * L
      } else if (type == "prec"){
        pattern <- matrix(modMatrices$kappa_c[[1]][idx, idx, g] != 0, p, p)
        K <- blockPrecision(Sk, pattern | t(pattern))
        modMatrices$kappa_c$start[idx, idx, g] <- pattern * K
      } else {
        pattern <- matrix(modMatrices$omega_c[[1]][idx, idx, g] != 0, p, p)
        K <- blockPrecision(Sk, pattern | t(pattern))
        omega <- -cov2cor(K)
        diag(omega) <- 0
        modMatrices$omega_c$start[idx, idx, g] <- (pattern | t(pattern)) * omega
        inModel <- diag(matrix(modMatrices$delta_c[[1]][idx, idx, g] != 0, p, p))
        modMatrices$delta_c$start[idx, idx, g] <- diag(inModel / sqrt(diag(K)), p, p)
      }
    }
  }
  modMatrices
}

# Composite variable model:
cvm <- function(
  data, # Dataset
  lambda, # Factor loadings (variables x constructs); only for the common factors
  weights, # Composite weights (variables x constructs); only for the composites
  latent = c("cov","ggm"), # Type of the (residual) covariances of the constructs
  residual = c("cov","chol","prec","ggm","cor"), # Type of the residual covariances of the factor indicators
  composite = c("cov","chol","prec","ggm","cor"), # Type of the covariances of the composite indicators (the "_c" matrices)

  # Latent matrices:
  sigma_zeta = "full",
  omega_zeta = "full", # Partial correlations
  delta_zeta = "full", # Scaling

  # Residual matrices (factor indicators):
  sigma_epsilon = "diag",
  kappa_epsilon = "diag",
  omega_epsilon = "zero",
  lowertri_epsilon = "diag",
  delta_epsilon = "diag",
  rho_epsilon = "zero",
  SD_epsilon = "diag",

  # Composite indicator matrices; "full" means full within every composite:
  sigma_c = "full",
  kappa_c = "full",
  omega_c = "full",
  lowertri_c = "full",
  delta_c = "full",
  rho_c = "full",
  SD_c = "full",

  # Beta:
  beta = "zero",

  # Mean structure:
  nu,
  nu_eta,

  # Identification:
  identify = TRUE,
  identification = c("loadings","variance"),
  bad_marker_crit = 1, # A composite is scaled by the indicator with the largest (standardized) starting weight if the starting weight of its first indicator is smaller than this proportion of the largest. 1: always the indicator with the largest starting weight; 0: always the first indicator
  sumscores = FALSE, # Fix all weights to 1, so that every composite is the sum score of its indicators?

  # Rest:
  vars, # character indicating the variables Extracted if missing from data - group variable
  latents, # Name of the constructs (latent variables and composites)
  groups, # deprecated, use groupvar instead
  groupvar, # grouping variable (character column name or vector of group names)
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "auto",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  estimator = "default",
  optimizer,
  storedata = FALSE,
  WLS.W,
  covtype = c("choose","ML","UB"),
  standardize = c("none","z","quantile"),
  sampleStats,
  verbose = FALSE,
  bootstrap = FALSE,
  boot_sub,
  boot_resample,
  # Penalized ML arguments (estimator = "PML" / "PFIML"):
  penalty_lambda = NA,  # Penalty strength (NA = auto-select via EBIC grid search)
  penalty_alpha = 1,   # Elastic net mixing: 1 = LASSO, 0 = ridge
  penalize_matrices  # Character vector of matrix names to penalize. Default: defaultPenalizeMatrices()
){
  # Experimental-feature note, gated at the shared cutoff (silent from 0.19);
  # fires once per session, unconditional of verbose (cnm() calls cvm()):
  experimentalWarning("cvm() / cnm() composite variable model family")

  # Standardize input arguments:
  si <- standardize_input(
    data = if(missing(data)) NULL else data,
    covs = if(missing(covs)) NULL else covs,
    cors = NULL,
    nobs = if(missing(nobs)) NULL else nobs,
    corinput = NULL,
    groups = if(missing(groups)) NULL else groups,
    groupvar = if(missing(groupvar)) NULL else groupvar,
    family = "cvm", caller = "cvm()", estimator = estimator
  )
  if (!is.null(si$data)) data <- si$data
  if (!is.null(si$covs)) covs <- si$covs
  if (!is.null(si$nobs)) nobs <- si$nobs
  if (!is.null(si$groups)) groups <- si$groups

  rawts <- FALSE

  # Estimator:
  if (estimator == "default"){
    estimator <- "ML"
  }
  if (!estimator %in% c("ML","FIML","ULS","WLS","DWLS","PML","PFIML")){
    stop("cvm() currently supports the estimators 'ML', 'FIML', 'ULS', 'WLS', 'DWLS', 'PML' and 'PFIML'.")
  }

  # Type:
  latent <- match.arg(latent)
  residual <- match.arg(residual)
  composite <- match.arg(composite)

  # Identification:
  identification <- match.arg(identification)

  # Sum scores (all weights fixed to 1):
  if (!is.logical(sumscores) || length(sumscores) != 1 || is.na(sumscores)){
    stop("'sumscores' must be TRUE or FALSE.")
  }

  # WLS weights:
  if (missing(WLS.W)){
    WLS.W <- ifelse(!estimator %in% c("WLS","ULS","DWLS"), "none",
                    switch(estimator,
                           "WLS" = "full",
                           "ULS" = "identity",
                           "DWLS" = "diag"
                    ))
  }

  # Auto-detect missing data handling:
  if (missing == "auto") {
    if (!missing(data)) {
      if (missing(vars)) {
        check_vars <- if (!is.null(colnames(data))) colnames(data) else seq_len(ncol(data))
        if (!missing(groups) && is.character(groups) && length(groups) == 1) check_vars <- setdiff(check_vars, groups)
      } else {
        check_vars <- vars
      }
      has_missing <- any(is.na(data[, check_vars, drop = FALSE]))
      if (has_missing) {
        if (estimator == "ML") {
          estimator <- "FIML"
        } else if (estimator == "PML") {
          estimator <- "PFIML"
        } else {
          missing <- "listwise"
        }
      } else {
        missing <- "listwise"
      }
    } else {
      missing <- "listwise"
    }
  }

  # Obtain sample stats:
  if (missing(sampleStats)){
    sampleStats <- samplestats(data = data,
                               vars = vars,
                               groups = groups,
                               covs = covs,
                               means = means,
                               nobs = nobs,
                               missing = ifelse(estimator %in% c("FIML","PFIML"),"pairwise",missing),
                               rawts = rawts,
                               fimldata = estimator %in% c("FIML","PFIML"),
                               storedata = storedata,
                               covtype = covtype,
                               weightsmatrix = WLS.W,
                               verbose = verbose,
                               standardize = standardize,
                               bootstrap = bootstrap,
                               boot_sub = boot_sub,
                               boot_resample = boot_resample)
  }

  # No correlation input and always a mean structure:
  sampleStats@corinput <- FALSE
  meanstructure <- TRUE

  # Number of variables:
  nNode <- nrow(sampleStats@variables)
  varnames <- sampleStats@variables$label

  # Generate model object:
  model <- generate_psychonetrics(model = "cvm",sample = sampleStats,computed = FALSE,
                                  equal = equal,identification=identification,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  rawts = rawts, types = list(latent = latent, residual = residual, composite = composite, likelihood = "normal"),
                                  meanstructure = meanstructure,
                                  verbose=verbose)

  # Submodel:
  model@submodel <- ifelse(composite %in% c("ggm","prec"), "cnm", "cvm")

  # Number of groups:
  nGroup <- nrow(model@sample@groups)

  # Number of means:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  # Add number of observations:
  model@sample@nobs <-
    nNode * (nNode-1) / 2 * nGroup + # Off-diagonal covariances per group
    nNode * nGroup + # Variances
    meanstructure * nMeans # Means per group

  ### Measurement structure ###
  if (missing(lambda) && missing(weights)){
    stop("At least one of 'lambda' (factor loadings) and 'weights' (composite weights) must be supplied.")
  }
  if (!missing(lambda) && is.character(lambda)) stop("'lambda' must be a numeric matrix.")
  if (!missing(weights) && is.character(weights)) stop("'weights' must be a numeric matrix.")
  if (missing(lambda)){
    lambda <- matrix(0, nNode, dim(weights)[2])
  }
  if (missing(weights)){
    weights <- matrix(0, nNode, dim(lambda)[2])
  }
  if (dim(lambda)[1] != nNode || dim(weights)[1] != nNode){
    stop("'lambda' and 'weights' must have one row per observed variable.")
  }
  if (dim(lambda)[2] != dim(weights)[2]){
    stop("'lambda' and 'weights' must have the same number of columns: one per construct (latent variable or composite).")
  }

  # Number of constructs:
  nLatent <- dim(lambda)[2]

  # Fix the pattern matrices:
  lambda <- fixMatrix(1 * lambda, nGroup = nGroup, equal = "lambda" %in% equal)
  weights <- fixMatrix(1 * weights, nGroup = nGroup, equal = "weights" %in% equal)

  # Structure (the same in every group):
  lamAny <- apply(lambda != 0, c(1,2), any)
  wAny <- apply(weights != 0, c(1,2), any)

  # Which constructs are composites?
  comp <- colSums(wAny) > 0
  fac <- colSums(lamAny) > 0
  if (any(comp & fac)){
    stop("A construct cannot have both factor loadings ('lambda') and composite weights ('weights'): a composite is formed by its indicators and has no reflective indicators. ",
         "To let an observed variable be predicted by (or predict) a composite, add that variable as a construct of its own with a single indicator and specify the effect in 'beta': ",
         "a latent variable with one loading in 'lambda' (its residual variance is then fixed to zero), or equivalently a composite with one weight in 'weights'.")
  }
  if (any(!comp & !fac)){
    stop("Every construct needs at least one indicator in 'lambda' or 'weights'.")
  }

  # To which composite does every variable belong (0: not a composite indicator)?
  if (any(rowSums(wAny) > 1)){
    stop("A variable can only be an indicator of a single composite.")
  }
  block <- as.vector(wAny %*% seq_len(nLatent))
  if (any(block != 0 & rowSums(lamAny) > 0)){
    stop("A variable cannot be both a factor indicator ('lambda') and a composite indicator ('weights').")
  }

  # Indices:
  compInd <- which(block != 0) # Composite indicators
  resInd <- which(block == 0) # Factor indicators (variables with a residual)
  nComp <- length(compInd)
  nRes <- length(resInd)
  compIdx <- which(comp)
  facIdx <- which(!comp)

  # Names of the constructs:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of constructs in the model.")
  }

  # The indicators of a composite may not be linearly dependent (error), and
  # a warning is given when they are nearly so:
  cvm_check_collinearity(model@sample@covs, block, compIdx, latents, sumscores = sumscores)

  ### Starting values ###
  lambdaStart <- vector("list", nGroup)
  weightsStart <- vector("list", nGroup)
  expLatSigma <- vector("list", nGroup)
  expResidSigma <- vector("list", nGroup)
  expCompSigma <- vector("list", nGroup)

  # Pattern of the structural model (needed to form the starting values):
  betaFixed <- fixMatrix(beta, nGroup = nGroup, nrows = nLatent, ncols = nLatent,
                         equal = "beta" %in% equal, diag0 = TRUE)

  # Which constructs are related to each other in the model?
  latentArg <- if (latent == "cov") sigma_zeta else omega_zeta
  if (is.character(latentArg)){
    covAny <- matrix(latentArg[1] == "full", nLatent, nLatent)
  } else {
    latentArr <- if (is.list(latentArg)) simplify2array(latentArg) else latentArg
    covAny <- if (length(dim(latentArr)) == 3) apply(latentArr != 0, c(1,2), any) else latentArr != 0
    covAny <- covAny | t(covAny)
  }
  betaAny <- apply(betaFixed != 0, c(1,2), any)
  adjacency <- betaAny | t(betaAny) | covAny
  diag(adjacency) <- FALSE

  # Is the construct part of the model saturated (all constructs are composites
  # that are all related to each other, without effects in beta)? The starting
  # weights can then be improved, see cvm_genvar_weights():
  saturatedConstructs <- nLatent > 1 && length(facIdx) == 0 && !any(betaAny) &&
    all(adjacency[upper.tri(adjacency)])

  # Factor part, using the lvm starting values on the factor indicators:
  if (length(facIdx) > 0 && nRes > 0){
    facSetup <- matrixsetup_lambda(lambda[resInd, facIdx, , drop = FALSE],
                                   expcov = lapply(model@sample@covs, function(S) as.matrix(S)[resInd, resInd, drop = FALSE]),
                                   nGroup = nGroup, equal = FALSE,
                                   observednames = varnames[resInd], latentnames = latents[facIdx],
                                   sampletable = sampleStats, identification = identification)
  }

  # First pass: starting values of the factor part and the (unit variance)
  # weights of the composites, per group:
  Sstart <- vector("list", nGroup)
  proxyStart <- vector("list", nGroup)
  facCovStart <- vector("list", nGroup)
  for (g in seq_len(nGroup)){
    S <- as.matrix(spectralshift(model@sample@covs[[g]]))
    Sstart[[g]] <- S

    lambdaStart[[g]] <- matrix(0, nNode, nLatent)
    facProxy <- matrix(0, nNode, nLatent) # weights that give a proxy score per factor

    if (length(facIdx) > 0 && nRes > 0){
      lambdaStart[[g]][resInd, facIdx] <- facSetup$start[,,g]
      expResidSigma[[g]] <- as.matrix(facSetup$sigma_epsilon_start[,,g])
      facCovStart[[g]] <- as.matrix(facSetup$sigma_zeta_start[,,g])
      # Least squares factor scores as proxies:
      facProxy[resInd, facIdx] <- t(MASS::ginv(as.matrix(facSetup$start[,,g])))
    } else if (nRes > 0){
      expResidSigma[[g]] <- diag(diag(S)[resInd], nRes, nRes)
    }

    if (sumscores){
      # Unit weights: every composite is the sum score of its indicators:
      proxyStart[[g]] <- facProxy
      for (k in compIdx) proxyStart[[g]][block == k, k] <- 1
    } else {
      proxyStart[[g]] <- cvm_start_weights(S, block, comp, facProxy, adjacency, genvar = saturatedConstructs)
    }

    # Covariances of the composite indicators:
    if (nComp > 0){
      expCompSigma[[g]] <- as.matrix(model@sample@covs[[g]])[compInd, compInd, drop = FALSE]
    }
  }

  # Which indicator scales (is the marker of) every composite? The first
  # indicator, unless its (standardized) starting weight is smaller than
  # 'bad_marker_crit' times the largest starting weight of the composite: the
  # indicator with the largest starting weight is then used instead. With the
  # default (1) a composite is always scaled by the indicator with the largest
  # starting weight. Scaling by an indicator with a small weight makes the other
  # weights large, and if that weight has another sign than the other weights
  # the composite is reversed; the optimizer then stopped before the optimum
  # in examples. Not used with sumscores = TRUE (all weights are then fixed).
  marker <- rep(NA_integer_, nLatent)
  rescaled <- character(0)
  for (k in compIdx){
    idx <- which(block == k)
    marker[k] <- idx[1]
    if (identify && !sumscores && length(idx) > 1 && is.numeric(bad_marker_crit) && bad_marker_crit > 0){
      # Standardized weights (weights times the indicator SD), averaged over groups:
      stdw <- rowMeans(do.call(cbind, lapply(seq_len(nGroup), function(g){
        abs(proxyStart[[g]][idx, k]) * sqrt(diag(Sstart[[g]])[idx])
      })))
      if (all(is.finite(stdw)) && max(stdw) > 0 && stdw[1] < bad_marker_crit * max(stdw)){
        marker[k] <- idx[which.max(stdw)]
        rescaled <- c(rescaled, paste0(latents[k], " (", varnames[marker[k]], ")"))
      }
    }
  }
  if (length(rescaled) > 0){
    message(paste0(
      if (bad_marker_crit >= 1){
        paste0("Composites are scaled by fixing the weight of the indicator with the largest (standardized) starting weight to 1. ",
               "For the following composite(s) this is not the first indicator (scaling indicator in parentheses): ")
      } else {
        paste0("The (standardized) starting weight of the first indicator of the following composite(s) is smaller than ", format(bad_marker_crit),
               " times the largest starting weight of the composite. These composites are scaled by fixing the weight of the indicator ",
               "with the largest starting weight to 1 instead (in parentheses): ")
      },
      paste(rescaled, collapse = ", "),
      ". Set 'bad_marker_crit = 0' to always use the first indicator."))
  }

  # Second pass: weights relative to the marker, and the covariances of the constructs:
  for (g in seq_len(nGroup)){
    S <- Sstart[[g]]
    sds <- sqrt(diag(S))
    proxy <- proxyStart[[g]]
    weightsStart[[g]] <- matrix(0, nNode, nLatent)

    for (k in compIdx){
      idx <- which(block == k)
      w <- proxy[idx, k] / proxy[marker[k], k]
      if (sumscores){
        w <- rep(1, length(idx))
      } else if (any(!is.finite(w)) || max(abs(w * sds[idx])) > 25 * sds[marker[k]]){
        # Equal standardized weights:
        w <- sds[marker[k]] / sds[idx]
      }
      weightsStart[[g]][idx, k] <- w
      proxy[idx, k] <- w
    }

    # Covariances of the constructs:
    latCov <- t(proxy) %*% S %*% proxy
    if (!is.null(facCovStart[[g]])){
      latCov[facIdx, facIdx] <- facCovStart[[g]]
    }
    expLatSigma[[g]] <- as.matrix(spectralshift(latCov))
  }

  ### Model matrices ###
  modMatrices <- list()

  # Intercepts of the factor indicators and means of the composite indicators:
  modMatrices$nu <- matrixsetup_mu(nu,nNode = nNode,nGroup = nGroup,labels = varnames,equal = "nu" %in% equal,
                                   expmeans = model@sample@means, sampletable = sampleStats, name = "nu")

  # Intercepts of the constructs:
  modMatrices$nu_eta <- matrixsetup_mu(nu_eta,nNode = nLatent,nGroup = nGroup,labels = latents,equal = "nu_eta" %in% equal,
                                       expmeans = lapply(seq_len(nGroup),function(x)rep(0,nLatent)), sampletable = sampleStats, name = "nu_eta")

  # Factor loadings:
  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov = model@sample@covs, nGroup = nGroup, equal = FALSE,
                                           observednames = varnames, latentnames = latents,
                                           sampletable = sampleStats, identification = identification,
                                           start = lambdaStart)

  # Composite weights:
  weightsStartArray <- weights
  for (g in seq_len(nGroup)){
    weightsStartArray[,,g] <- (weights[,,g] != 0) * weightsStart[[g]]
  }
  modMatrices$weights <- list(weights,
                              mat = "weights",
                              op = "~>",
                              rownames = varnames,
                              colnames = latents,
                              sparse = TRUE,
                              start = weightsStartArray,
                              sampletable = sampleStats)

  # Beta, with OLS based starting values:
  beta_start_list <- NULL
  if (any(betaFixed != 0)){
    beta_start_list <- cvm_beta_start(betaFixed, expLatSigma)
  }
  if (!is.null(beta_start_list)){
    modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                         labels = latents, sampletable = sampleStats,
                                         equal = "beta" %in% equal,
                                         start = beta_start_list)
    # From total to residual covariances:
    for (g in seq_len(nGroup)){
      IminB <- diag(nLatent) - beta_start_list[[g]]
      expLatSigma[[g]] <- as.matrix(spectralshift(IminB %*% expLatSigma[[g]] %*% t(IminB)))
    }
  } else {
    modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                         labels = latents, sampletable = sampleStats,
                                         equal = "beta" %in% equal)
  }

  # Latent varcov:
  if (latent == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta,
                                                name = "sigma_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = "sigma_zeta" %in% equal, sampletable = sampleStats,
                                                beta = modMatrices$beta[[1]])
  } else if (latent == "ggm"){
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta,
                                                name = "omega_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = "omega_zeta" %in% equal, sampletable = sampleStats,
                                                beta = modMatrices$beta[[1]])

    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta,
                                                name = "delta_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = "delta_zeta" %in% equal, sampletable = sampleStats,
                                                omegaStart =  modMatrices$omega_zeta$start)
  }

  # Helper to add a typed covariance block:
  addBlock <- function(modMatrices, type, suffix, args, expcov, n, labels, ids, pattern = identity, kinds){
    setup <- function(fun, what, kind, ...){
      name <- paste0(what, "_", suffix)
      res <- fun(pattern(args[[what]], name, kind), name = name, expcov = expcov, nNode = n, nGroup = nGroup,
                 labels = labels, equal = name %in% equal, sampletable = sampleStats, ...)
      res$rowid <- ids
      res$colid <- ids
      res
    }
    if (type == "cov"){
      modMatrices[[paste0("sigma_", suffix)]] <- setup(matrixsetup_sigma, "sigma", "symmetric")
    } else if (type == "chol"){
      modMatrices[[paste0("lowertri_", suffix)]] <- setup(matrixsetup_lowertri, "lowertri", "symmetric")
    } else if (type == "ggm"){
      modMatrices[[paste0("omega_", suffix)]] <- setup(matrixsetup_omega, "omega", "offdiagonal")
      modMatrices[[paste0("delta_", suffix)]] <- setup(matrixsetup_delta, "delta", "diagonal",
                                                         omegaStart = modMatrices[[paste0("omega_", suffix)]]$start)
    } else if (type == "prec"){
      modMatrices[[paste0("kappa_", suffix)]] <- setup(matrixsetup_kappa, "kappa", "symmetric")
    } else if (type == "cor"){
      modMatrices[[paste0("rho_", suffix)]] <- setup(matrixsetup_rho, "rho", "offdiagonal")
      modMatrices[[paste0("SD_", suffix)]] <- setup(matrixsetup_SD, "SD", "diagonal")
    }
    modMatrices
  }

  # Residual varcov of the factor indicators:
  if (nRes > 0){
    modMatrices <- addBlock(modMatrices, residual, "epsilon",
                            args = list(sigma = sigma_epsilon, lowertri = lowertri_epsilon, omega = omega_epsilon,
                                        delta = delta_epsilon, kappa = kappa_epsilon, rho = rho_epsilon, SD = SD_epsilon),
                            expcov = expResidSigma, n = nRes, labels = varnames[resInd], ids = resInd,
                            pattern = function(x, name, kind) x)
  }

  # Varcov of the composite indicators (block diagonal over composites):
  if (nComp > 0){
    blockComp <- block[compInd]
    # The starting values of the network, precision and Cholesky types are
    # formed per composite below (cvm_block_start), so the shared setup
    # functions only need the variances here:
    expCompSetup <- expCompSigma
    if (composite %in% c("ggm", "prec", "chol")){
      expCompSetup <- lapply(expCompSigma, function(S){
        if (any(!is.finite(S))) S else diag(diag(S), nrow(S), ncol(S))
      })
    }
    modMatrices <- addBlock(modMatrices, composite, "c",
                            args = list(sigma = sigma_c, lowertri = lowertri_c, omega = omega_c,
                                        delta = delta_c, kappa = kappa_c, rho = rho_c, SD = SD_c),
                            expcov = expCompSetup, n = nComp, labels = varnames[compInd], ids = compInd,
                            pattern = function(x, name, kind) cvm_blockpattern(x, blockComp, nGroup, name, kind))
    modMatrices <- cvm_block_start(modMatrices, composite, expCompSigma, blockComp, nGroup)
  }

  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)

  # Composites that are scaled by another indicator than the first: fix the
  # weight of that indicator to 1 (the automatic identification then leaves the
  # composite as is):
  for (k in compIdx){
    if (!is.na(marker[k]) && marker[k] != which(block == k)[1]){
      rows <- which(pars$partable$matrix == "weights" & pars$partable$row == marker[k] & pars$partable$col == k)
      pars$partable$est[rows] <- 1
      pars$partable$par[rows] <- 0
      pars$partable$fixed[rows] <- TRUE
      pars$partable$identified[rows] <- TRUE
    }
  }

  # Sum scores: fix all weights to 1. These are restrictions of the model
  # rather than identification constraints (as when weights are fixed with
  # fixpar()), so they are not marked as 'identified' and modification indices
  # are computed for them:
  if (sumscores){
    rows <- which(pars$partable$matrix == "weights" & !pars$partable$fixed)
    pars$partable$est[rows] <- 1
    pars$partable$par[rows] <- 0
    pars$partable$fixed[rows] <- TRUE
  }
  pars$partable <- parRelabel(pars$partable)

  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable

  # Extra matrices:
  Lmat <- psychonetrics::eliminationMatrix(nNode)
  extramatrices <- list(
    D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    Deta = psychonetrics::duplicationMatrix(nLatent), # non-strict duplciation matrix
    L = Lmat, # Elinimation matrix
    L_eta = psychonetrics::eliminationMatrix(nLatent), # Elinimation matrix
    Dstar_eta = psychonetrics::duplicationMatrix(nLatent,diag = FALSE), # Strict duplicaton matrix
    In = as(diag(nNode),"dMatrix"), # Identity of dim n
    Inlatent = as(diag(nLatent),"dMatrix"),
    C = commutationMatrix(nNode, nLatent),
    Cbeta = commutationMatrix(nLatent, nLatent),
    Aeta = psychonetrics::diagonalizationMatrix(nLatent),
    compvec = matrix(as.numeric(comp), ncol = 1), # Which constructs are composites?
    blockvec = matrix(block, ncol = 1) # Composite of every variable (0: not a composite indicator)
  )

  # Selection matrix and algebra helpers of a block of variables:
  blockMatrices <- function(ind, suffix){
    n <- length(ind)
    E <- Matrix::sparseMatrix(i = ind, j = seq_len(n), x = 1, dims = c(nNode, n))
    E <- as(E, "dMatrix")
    res <- list(
      E = E,
      L = psychonetrics::eliminationMatrix(n),
      D = psychonetrics::duplicationMatrix(n),
      Dstar = psychonetrics::duplicationMatrix(n, diag = FALSE),
      In = as(diag(n),"dMatrix"),
      C_chol = commutationMatrix(n, n),
      A = psychonetrics::diagonalizationMatrix(n)
    )
    names(res) <- paste0(names(res), "_", suffix)
    res
  }
  if (nRes > 0){
    extramatrices <- c(extramatrices, blockMatrices(resInd, "res"))
    # Maps vech of the residual block to vech of the full matrix:
    extramatrices$Sel_res <- as(Lmat %*% kronecker(extramatrices$E_res, extramatrices$E_res) %*% extramatrices$D_res, "dMatrix")
  }
  if (nComp > 0){
    extramatrices <- c(extramatrices, blockMatrices(compInd, "comp"))
    # Maps vec of the composite indicator block to vech of the full matrix:
    extramatrices$LEE_comp <- as(Lmat %*% kronecker(extramatrices$E_comp, extramatrices$E_comp), "dMatrix")
  }
  model@extramatrices <- extramatrices

  ### Baseline model ###
  # The baseline and saturated models always use the unpenalized estimator:
  base_est <- switch(estimator, "PML" = "ML", "PFIML" = "FIML", estimator)
  if (baseline_saturated){
    model@baseline_saturated$baseline <- varcov(data,
                                                type = "chol",
                                                lowertri = "diag",
                                                vars = vars,
                                                groups = groups,
                                                covs = covs,
                                                means = means,
                                                nobs = nobs,
                                                missing = missing,
                                                equal = equal,
                                                estimator = base_est,
                                                baseline_saturated = FALSE, sampleStats = sampleStats)

    ### Saturated model ###
    model@baseline_saturated$saturated <- varcov(data,
                                                 type = "chol",
                                                 lowertri = "full",
                                                 vars = vars,
                                                 groups = groups,
                                                 covs = covs,
                                                 means = means,
                                                 nobs = nobs,
                                                 missing = missing,
                                                 estimator = base_est,
                                                 baseline_saturated = FALSE, sampleStats = sampleStats)

    # if not FIML, Treat as computed:
    if (base_est != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
    }
  }

  # Identify model (also marks the computed entries and structural zeros, so
  # this is always needed):
  model <- identify_cvm(model, identify = identify)

  # Make sure the starting values are admissible (positive computed variances
  # of composites and positive definite implied covariances):
  model <- cvm_admissible_start(model)

  # Start the intercepts of the factor indicators such that the implied means
  # equal the sample means, given the implied means of the constructs at the
  # starting values (see cvm_start_intercepts):
  model <- cvm_start_intercepts(model)

  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)

  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # Setup PML/PFIML penalization (by default of the network structures in the
  # model, see defaultPenalizeMatrices):
  if (estimator %in% c("PML", "PFIML")) {
    model@penalty <- list(lambda = penalty_lambda, alpha = penalty_alpha)
    pen_mats <- if (missing(penalize_matrices)) defaultPenalizeMatrices(model) else penalize_matrices
    if (length(pen_mats) == 0){
      warning("No matrices to penalize: by default only network structures (omega_c / kappa_c, omega_epsilon / kappa_epsilon, omega_zeta) are penalized in cvm(). ",
              "Use composite = 'ggm' (or cnm()), residual = 'ggm' or latent = 'ggm', or set 'penalize_matrices'.")
    } else {
      model <- penalize(model, matrix = pen_mats, lambda = penalty_lambda, log = FALSE)
    }
  }

  # Return model:
  return(model)
}

# Composite network model (cnm): a cvm in which the covariances of the
# composite indicators are modeled as a Gaussian graphical model, resulting in
# one network of indicators per composite (omega_c). By default the constructs
# are modeled as a network as well (latent = "ggm", omega_zeta).
cnm <- function(...){
  # Which arguments of cvm() are supplied (also when abbreviated or unnamed)?
  given <- tryCatch(names(as.list(match.call(cvm, sys.call()))[-1]),
                    error = function(e) names(list(...)))
  if ("composite" %in% given){
    stop("'composite' cannot be set in cnm(): the composite indicators are modeled as a network (composite = 'ggm'). Use cvm() for other types.")
  }
  if ("latent" %in% given){
    return(cvm(..., composite = "ggm"))
  }
  if ("sigma_zeta" %in% given){
    stop("'sigma_zeta' is not used by default in cnm(): the constructs are modeled as a network (latent = 'ggm'), which is specified with 'omega_zeta'. Set latent = 'cov' to model the covariances of the constructs with 'sigma_zeta'.")
  }
  cvm(..., latent = "ggm", composite = "ggm")
}
