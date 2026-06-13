# Convert a lavaan model into an equivalent psychonetrics lvm.
#
# Only the lavaan public API is used (lavInspect, parTable, lavNames, sem); no
# lavaan source code is ported. See man/fromlavaan.Rd for the conventions that
# are guaranteed to be reproduced and for the list of unsupported features.
fromlavaan <- function(x,
                       data = NULL,
                       run = FALSE,
                       estimator = "default",
                       baseline_saturated = TRUE,
                       verbose = FALSE,
                       ...){
  if (verbose) experimentalWarning("fromlavaan()")

  # Obtain a (possibly unfitted) lavaan object:
  if (is(x, "lavaan")){
    fit <- x
  } else {
    # Treat as model syntax; data is required:
    if (is.null(data)){
      stop("When 'x' is lavaan model syntax, 'data' must be supplied.")
    }
    fit <- lavaan::sem(model = x, data = data, meanstructure = TRUE,
                       fixed.x = FALSE, do.fit = FALSE)
  }

  pt  <- lavaan::parTable(fit)
  opt <- lavaan::lavInspect(fit, "options")

  ## ---- Scope guards (validated in the prototypes) -------------------
  if (lavaan::lavInspect(fit, "nlevels") > 1){
    stop("Two-level (cluster=/level:) lavaan models cannot be converted to psychonetrics.")
  }
  if (any(pt$op == ":=")){
    stop("Defined parameters (':=') cannot be converted to psychonetrics.")
  }
  if (any(pt$op %in% c("<", ">"))){
    stop("Inequality constraints ('<', '>') cannot be converted to psychonetrics.")
  }
  if (any(pt$op == "<~")){
    stop("Formative indicators ('<~') cannot be converted to psychonetrics.")
  }
  if (isTRUE(opt$conditional.x)){
    stop("Models with conditional.x = TRUE cannot be converted to psychonetrics.")
  }
  if (length(lavaan::lavNames(fit, "ov.ord")) > 0){
    stop("Categorical (ordinal) endogenous variables are not yet supported by fromlavaan().")
  }

  # Estimator detection / FIML:
  if (isTRUE(opt$missing %in% c("ml", "fiml", "ml.x"))){
    if (!identical(estimator, "default") && !identical(estimator, "FIML")){
      warning("Missing data detected; overriding estimator to 'FIML'.")
    }
    estimator <- "FIML"
  } else if (identical(estimator, "default")){
    estimator <- "ML"
  }

  # Likelihood / robustness warnings:
  if (identical(opt$likelihood, "wishart")){
    warning("lavaan used likelihood = 'wishart' (n - 1 normalization); psychonetrics ",
            "uses the normal (n) likelihood. Estimates match but the log-likelihood ",
            "and any (n - 1)-based fit measures will differ.")
  }
  # Robust SEs/tests are not transferred. "none" (e.g. when do.fit = FALSE) and
  # "standard" are not robust, so they should not trigger the warning.
  robust_se   <- !is.null(opt$se)   && !all(opt$se %in% c("standard", "none"))
  robust_test <- !is.null(opt$test) && !all(opt$test %in% c("standard", "none"))
  if (robust_se || robust_test){
    warning("Robust standard errors / test statistics are not transferred; the model ",
            "is converted as plain ML.")
  }

  nG            <- lavaan::lavInspect(fit, "ngroups")
  meanstructure <- isTRUE(lavaan::lavInspect(fit, "meanstructure"))

  # fixed.x = TRUE with exogenous covariates: loud warning.
  if (isTRUE(opt$fixed.x) && length(lavaan::lavNames(fit, "ov.x")) > 0){
    warning("lavaan used fixed.x = TRUE with exogenous covariates. psychonetrics treats ",
            "these as endogenous (random) variables: parameter estimates, standard errors ",
            "and the chi-square will match lavaan, but df-based fit measures and the ",
            "log-likelihood will NOT (lavaan conditions on the exogenous covariates).")
  }

  ## ---- Free / est matrices (lavaan LISREL representation) -----------
  FREE <- lavaan::lavInspect(fit, "free")
  EST  <- lavaan::lavInspect(fit, "est")
  if (nG == 1){ FREE <- list(FREE); EST <- list(EST) }

  ov  <- rownames(FREE[[1]]$lambda)
  eta <- colnames(FREE[[1]]$lambda)
  # Avoid latent/observed name collisions (observed structural variables enter
  # lambda as phantom latent columns):
  eta_pn <- ifelse(eta %in% ov, paste0("eta_", eta), eta)
  hasBeta <- !is.null(FREE[[1]]$beta)

  ## ---- Equality classes from the parameter table --------------------
  # Union-find over lavaan free indices 1..max(pt$free).
  freerows <- pt[pt$free > 0, ]
  maxfree  <- max(pt$free)
  cl <- seq_len(maxfree)
  findroot <- function(a){ while (cl[a] != a) a <- cl[a]; a }
  uni <- function(a, b){ ra <- findroot(a); rb <- findroot(b); if (ra != rb) cl[rb] <<- ra; invisible(NULL) }
  # (a) shared (nonempty) labels:
  labs <- freerows$label[freerows$label != ""]
  for (L in unique(labs)){
    f <- freerows$free[freerows$label == L]
    if (length(f) > 1) for (k in f[-1]) uni(f[1], k)
  }
  # (b) simple "==" rows: resolve each side against plabel, then label.
  eqs <- pt[pt$op == "==", ]
  resolve <- function(s){
    f <- freerows$free[freerows$plabel == s]
    if (length(f) == 0) f <- freerows$free[freerows$label == s]
    f
  }
  if (nrow(eqs) > 0) for (i in seq_len(nrow(eqs))){
    fl <- resolve(eqs$lhs[i]); fr <- resolve(eqs$rhs[i])
    if (length(fl) == 0 || length(fr) == 0){
      stop("General (non-simple) equality constraints cannot be converted to psychonetrics: ",
           eqs$lhs[i], " == ", eqs$rhs[i])
    }
    for (k in c(fl[-1], fr)) uni(fl[1], k)
  }
  # Resolve every index to its root so class labels are stable:
  cl <- vapply(cl, findroot, numeric(1))

  ## ---- Skeletons for lvm() ------------------------------------------
  skel <- function(mat){ # 1 where free OR fixed-at-nonzero, as a per-group 3-D array
    arr <- array(0, c(nrow(FREE[[1]][[mat]]), ncol(FREE[[1]][[mat]]), nG))
    for (g in 1:nG) arr[, , g] <- 1 * (FREE[[g]][[mat]] != 0 | EST[[g]][[mat]] != 0)
    arr
  }
  lambda_skel <- skel("lambda")
  psi_skel    <- skel("psi")
  theta_skel  <- skel("theta")
  beta_skel   <- if (hasBeta) skel("beta") else "zero"

  ## ---- Sample statistics --------------------------------------------
  ss <- lavaan::lavInspect(fit, "sampstat")
  if (nG == 1) ss <- list(ss)
  # Strip the "lavaan.matrix.symmetric"/"lavaan.vector" classes:
  covs <- lapply(ss, function(s){
    n <- nrow(s$cov)
    matrix(as.numeric(s$cov), n, n, dimnames = dimnames(s$cov))
  })
  means <- if (meanstructure) lapply(ss, function(s) as.numeric(s$mean)) else NULL
  nobs  <- unlist(lavaan::lavInspect(fit, "nobs"))

  args <- list(lambda = lambda_skel,
               latent = "cov", residual = "cov",
               sigma_zeta = psi_skel, sigma_epsilon = theta_skel,
               beta = beta_skel,
               vars = ov, latents = eta_pn,
               covs = covs, nobs = nobs, covtype = "ML",
               identify = FALSE,
               identification = if (isTRUE(opt$std.lv)) "variance" else "loadings",
               estimator = estimator,
               baseline_saturated = baseline_saturated,
               verbose = verbose)
  if (meanstructure) args$means <- means
  if (identical(estimator, "FIML")){
    # FIML needs raw data:
    dat <- lavaan::lavInspect(fit, "data")
    if (nG > 1){
      glabels <- lavaan::lavInspect(fit, "group.label")
      dat <- do.call(rbind, lapply(1:nG, function(g)
        data.frame(as.data.frame(dat[[g]]), .group = glabels[g],
                   stringsAsFactors = FALSE, check.names = FALSE)))
      args$groupvar <- ".group"
    } else {
      dat <- as.data.frame(dat)
    }
    args$data  <- dat
    args$covs  <- NULL
    args$means <- NULL
    args$nobs  <- NULL
  }
  # Pass through user arguments to lvm():
  args <- modifyList(args, list(...))
  mod <- do.call(lvm, args)

  ## ---- Parameter-table surgery --------------------------------------
  P <- mod@parameters
  matmap <- c(lambda = "lambda", theta = "sigma_epsilon", psi = "sigma_zeta",
              beta = "beta", nu = "nu", alpha = "nu_eta")
  sym <- c(theta = TRUE, psi = TRUE, lambda = FALSE, beta = FALSE,
           nu = FALSE, alpha = FALSE)
  # When lavaan has no meanstructure, leave nu saturated (excluded from the
  # reset and from the surgery): psychonetrics adds a saturated mean structure.
  lavmats <- if (meanstructure) names(matmap) else setdiff(names(matmap), "nu")
  # Pre-reset every managed matrix to fixed-at-zero; surgery re-frees below.
  rs <- P$matrix %in% unname(matmap[lavmats])
  P$par[rs] <- 0; P$fixed[rs] <- TRUE; P$est[rs] <- 0; P$identified[rs] <- FALSE

  for (g in 1:nG){
    for (m in lavmats){
      Fm <- FREE[[g]][[m]]; Em <- EST[[g]][[m]]
      if (is.null(Fm)) next
      pm <- matmap[[m]]
      for (i in seq_len(nrow(Fm))) for (j in seq_len(ncol(Fm))){
        if (sym[[m]] && j > i) next
        f <- Fm[i, j]; v <- Em[i, j]
        if (f == 0 && v == 0) next      # fixed at zero: already so
        idx <- which(P$matrix == pm & P$row == i & P$col == j & P$group_id == g)
        if (length(idx) != 1) stop("Internal error: no unique parameter row for ",
                                   pm, "[", i, ",", j, "] in group ", g, ".")
        P$est[idx] <- v
        if (f == 0){                     # fixed nonzero (marker, fixed.x, growth)
          P$fixed[idx] <- TRUE; P$par[idx] <- 0
          P$identified[idx] <- (pm == "lambda" && v == 1) ||
                               (pm == "sigma_zeta" && i == j && v == 1)
        } else {                         # free; par = equality class
          P$fixed[idx] <- FALSE; P$identified[idx] <- FALSE
          P$par[idx]   <- cl[f]
        }
      }
    }
  }

  if (!meanstructure){
    warning("lavaan model has no meanstructure; psychonetrics adds a saturated mean ",
            "structure (free nu, nu_eta fixed to 0). Chi-square, df, estimates and ",
            "standard errors match lavaan; the log-likelihood differs by a constant.")
  }

  P <- parRelabel(P)
  mod@parameters    <- P
  mod@modelmatrices <- formModelMatrices(mod)

  # Log:
  src <- if (is(x, "lavaan")) "lavaan object" else "lavaan model syntax"
  mod <- addLog(mod, paste0("Model created from a ", src, " via fromlavaan()."))
  mod@log[[1]] <- NULL # Remove the auto "Model created" entry.

  if (run) mod <- runmodel(mod)
  mod
}
