# Convert a lavaan model into an equivalent psychonetrics lvm, or, when the
# lavaan model contains composites ('<~'), an equivalent psychonetrics cvm.
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
  # Composites ('<~') are converted to the composite variable model, cvm():
  hasComposites <- any(pt$op == "<~")
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

  ## ---- Composites: structure and scope ------------------------------
  # compInd / resInd: observed variables that are / are not composite
  # indicators; blockOf: the composite (column) every variable belongs to.
  compInd <- integer(0); resInd <- seq_along(ov); blockOf <- rep(0, length(ov))
  autoFixed <- NULL
  if (hasComposites){
    if (is.null(FREE[[1]]$wmat)){
      stop("Formative indicators ('<~') can only be converted when lavaan models them as composites ",
           "(lavaan >= 0.7 with composites = TRUE).")
    }
    wAny <- Reduce(`|`, lapply(seq_len(nG), function(g) FREE[[g]]$wmat != 0 | EST[[g]]$wmat != 0))
    if (any(rowSums(wAny) > 1)){
      stop("A variable that is an indicator of more than one composite cannot be converted to psychonetrics.")
    }
    blockOf <- as.vector(wAny %*% seq_len(ncol(wAny)))
    compInd <- which(blockOf != 0)
    resInd  <- which(blockOf == 0)
    sameBlock <- outer(blockOf, blockOf, "==") & outer(blockOf != 0, blockOf != 0, "&")
    for (g in seq_len(nG)){
      thAny <- FREE[[g]]$theta != 0 | EST[[g]]$theta != 0
      thAny <- thAny | t(thAny)
      if (any(thAny[compInd, resInd])){
        stop("(Residual) covariances between an indicator of a composite and another observed variable ",
             "cannot be converted to psychonetrics.")
      }
      if (any(thAny[compInd, compInd, drop = FALSE] & !sameBlock[compInd, compInd, drop = FALSE])){
        stop("Covariances between indicators of different composites cannot be converted to psychonetrics: ",
             "in cvm() these are implied by the relations between the composites.")
      }
      # psychonetrics scales every composite with a fixed (marker) weight:
      for (k in which(colSums(wAny) > 0)){
        if (!any(FREE[[g]]$wmat[, k] == 0 & EST[[g]]$wmat[, k] != 0)){
          stop("Composite '", eta[k], "' has no fixed weight. psychonetrics scales a composite by fixing the weight ",
               "of one indicator; fix one of the weights in the lavaan model.")
        }
      }
    }
    # lavaan fixes the (co)variances of composite indicators to their sample
    # values by default (composites.cov = "fixed": free = 0 without a user
    # value). psychonetrics estimates them, which gives the same estimates
    # with complete data as these blocks are saturated:
    autoFixed <- lapply(seq_len(nG), function(g) matrix(FALSE, length(ov), length(ov)))
    rows <- which(pt$op == "~~" & pt$free == 0 & is.na(pt$ustart) &
                    pt$lhs %in% ov[compInd] & pt$rhs %in% ov[compInd])
    for (r in rows){
      g <- pt$group[r]; i <- match(pt$lhs[r], ov); j <- match(pt$rhs[r], ov)
      if (sameBlock[i, j]){
        autoFixed[[g]][i, j] <- autoFixed[[g]][j, i] <- TRUE
      }
    }
    if (any(unlist(autoFixed))){
      warning("lavaan fixed the (co)variances of the composite indicators to their sample values ",
              "(composites.cov = 'fixed', the lavaan default); psychonetrics estimates them. With complete data ",
              "the chi-square, degrees of freedom, estimates and the standard errors of weights, loadings and ",
              "regressions match lavaan; the number of parameters (AIC, BIC) and the standard errors of ",
              "(co)variances between composites do not. With missing data all results differ slightly. ",
              "Use composites.cov = 'free' in lavaan for an exact match.")
    }
  }
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
  if (hasComposites){
    weights_skel <- skel("wmat")
  }

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
  # Composites: a cvm with the residuals of the factor indicators in
  # sigma_epsilon and the (co)variances of the composite indicators in sigma_c:
  if (hasComposites){
    args$weights <- weights_skel
    args$composite <- "cov"
    args$sigma_c <- theta_skel[compInd, compInd, , drop = FALSE]
    if (length(resInd) > 0){
      args$sigma_epsilon <- theta_skel[resInd, resInd, , drop = FALSE]
    } else {
      args$sigma_epsilon <- NULL
    }
  }

  # Pass through user arguments to lvm() / cvm():
  args <- modifyList(args, list(...))
  mod <- do.call(if (hasComposites) cvm else lvm, args)

  ## ---- Parameter-table surgery --------------------------------------
  P <- mod@parameters
  matmap <- c(lambda = "lambda", theta = "sigma_epsilon", psi = "sigma_zeta",
              beta = "beta", nu = "nu", alpha = "nu_eta")
  if (hasComposites) matmap <- c(matmap, wmat = "weights")
  sym <- c(theta = TRUE, psi = TRUE, lambda = FALSE, beta = FALSE,
           nu = FALSE, alpha = FALSE, wmat = FALSE)
  # When lavaan has no meanstructure, leave nu saturated (excluded from the
  # reset and from the surgery): psychonetrics adds a saturated mean structure.
  lavmats <- if (meanstructure) names(matmap) else setdiff(names(matmap), "nu")
  # Pre-reset every managed matrix to fixed-at-zero; surgery re-frees below.
  managed <- unname(matmap[lavmats])
  if (hasComposites) managed <- c(managed, "sigma_c")
  rs <- P$matrix %in% managed
  est0 <- P$est   # starting values formed by lvm() / cvm()
  P$par[rs] <- 0; P$fixed[rs] <- TRUE; P$est[rs] <- 0; P$identified[rs] <- FALSE

  # Values of the free parameters: the lavaan estimates (or starting values).
  # Composite models are the exception when the lavaan model was not fitted or
  # did not converge: the lavaan values can then be far from a proper solution,
  # and the starting values of cvm() are used instead.
  useLavStart <- TRUE
  if (hasComposites){
    useLavStart <- isTRUE(lavaan::lavInspect(fit, "converged"))
    if (!useLavStart && isTRUE(opt$do.fit) && is(x, "lavaan")){
      warning("The lavaan model did not converge; its estimates are not used as starting values.")
    }
  }

  # Offset for the reassigned free-parameter indices. When nu is left saturated
  # (no meanstructure), its free rows KEEP the contiguous par indices that lvm()
  # assigned them (1..nNode); the surgery below re-frees the managed matrices
  # with the lavaan free indices cl[f], which also start at 1. Without an offset
  # the two index ranges overlap and parRelabel() (which merges rows sharing a
  # par value) would silently equate observed intercepts with loadings/variances
  # -- corrupting every single-group meanstructure = FALSE model. Shift the
  # reassigned indices past the largest retained (nu) par so the ranges are
  # disjoint; parRelabel() then compresses them back to a contiguous 1..K.
  par_offset <- if (meanstructure) 0L else max(c(0L, P$par[P$matrix == "nu"]))

  # Composites: entries that are computed in cvm() and never estimated (the
  # variances and intercepts of composites), and a counter for the composite
  # indicator (co)variances that lavaan fixed to their sample values but that
  # psychonetrics estimates:
  isComp <- rep(FALSE, length(eta))
  if (hasComposites) isComp <- colSums(wAny) > 0
  nextpar <- maxfree + par_offset

  for (g in 1:nG){
    for (m in lavmats){
      Fm <- FREE[[g]][[m]]; Em <- EST[[g]][[m]]
      if (is.null(Fm)) next
      for (i in seq_len(nrow(Fm))) for (j in seq_len(ncol(Fm))){
        if (sym[[m]] && j > i) next
        f <- Fm[i, j]; v <- Em[i, j]
        if (f == 0 && v == 0) next      # fixed at zero: already so

        # psychonetrics matrix and position of this element:
        pm <- matmap[[m]]; pr <- i; pc <- j
        computed <- FALSE; estimate <- FALSE
        if (hasComposites){
          if (m == "theta"){
            if (blockOf[i] != 0){       # composite indicators: sigma_c
              pm <- "sigma_c"; pr <- match(i, compInd); pc <- match(j, compInd)
              estimate <- f == 0 && autoFixed[[g]][i, j]
            } else {                    # other variables: sigma_epsilon
              pr <- match(i, resInd); pc <- match(j, resInd)
            }
          } else if ((m == "psi" && i == j && isComp[i]) || (m == "alpha" && isComp[i])){
            computed <- TRUE            # variance / intercept of a composite
          }
        }

        idx <- which(P$matrix == pm & P$row == pr & P$col == pc & P$group_id == g)
        if (length(idx) != 1) stop("Internal error: no unique parameter row for ",
                                   pm, "[", pr, ",", pc, "] in group ", g, ".")
        P$est[idx] <- if (useLavStart || (f == 0 && !estimate)) v else est0[idx]
        if (computed){                   # computed from the other parameters
          P$fixed[idx] <- TRUE; P$par[idx] <- 0; P$identified[idx] <- TRUE
        } else if (estimate){            # fixed to the sample value by lavaan
          nextpar <- nextpar + 1L
          P$fixed[idx] <- FALSE; P$identified[idx] <- FALSE
          P$par[idx]   <- nextpar
        } else if (f == 0){              # fixed nonzero (marker, fixed.x, growth)
          P$fixed[idx] <- TRUE; P$par[idx] <- 0
          P$identified[idx] <- (pm %in% c("lambda", "weights") && v == 1) ||
                               (pm == "sigma_zeta" && i == j && v == 1)
        } else {                         # free; par = equality class
          P$fixed[idx] <- FALSE; P$identified[idx] <- FALSE
          P$par[idx]   <- cl[f] + par_offset
        }
      }
    }
  }

  # Composites: the structural zeros and computed entries of the cvm framework
  # are never estimated, and are marked as such (as identify_cvm() does):
  if (hasComposites){
    struc <- cvm_structural_rows(mod, P)
    lock <- struc$zero | struc$intercept | struc$variance
    P$fixed[lock] <- TRUE; P$par[lock] <- 0; P$identified[lock] <- TRUE
    P$est[struc$zero] <- 0
    novalue <- struc$variance & (!is.finite(P$est) | P$est <= 0)
    P$est[novalue] <- 1

    # Starting values of cvm() are in the metric in which the first weight of a
    # composite is 1. If lavaan scales a composite differently (another marker
    # or a fixed value other than 1), bring the starting values of that
    # composite to the lavaan metric:
    if (!useLavStart){
      for (g in 1:nG) for (k in which(isComp)){
        wr <- which(P$matrix == "weights" & P$col == k & P$group_id == g & blockOf[P$row] == k)
        fx <- wr[P$fixed[wr] & P$est[wr] != 0]
        if (length(fx) == 0) next
        cc <- P$est[fx[1]] / est0[fx[1]]
        if (!is.finite(cc) || abs(cc) < 1e-3 || abs(cc) > 1e3 || abs(cc - 1) < 1e-12) next
        fr <- wr[!P$fixed[wr]]
        P$est[fr] <- P$est[fr] * cc
        isFreeG <- !P$fixed & P$group_id == g
        b_out <- which(isFreeG & P$matrix == "beta" & P$col == k & P$row != k)
        b_in  <- which(isFreeG & P$matrix == "beta" & P$row == k & P$col != k)
        P$est[b_out] <- P$est[b_out] / cc
        P$est[b_in]  <- P$est[b_in] * cc
        cv <- which(isFreeG & P$matrix == "sigma_zeta" & P$row != P$col & (P$row == k | P$col == k))
        P$est[cv] <- P$est[cv] * cc
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
  if (hasComposites){
    # Make sure the starting values are admissible (they are when taken from a
    # fitted lavaan model):
    mod <- cvm_admissible_start(mod)
    # With the psychonetrics starting values (lavaan model not fitted or not
    # converged), also start the intercepts of the factor indicators consistently
    # with the implied means of the constructs:
    if (!useLavStart) mod <- cvm_start_intercepts(mod)
  }
  mod@modelmatrices <- formModelMatrices(mod)

  # Log:
  src <- if (is(x, "lavaan")) "lavaan object" else "lavaan model syntax"
  mod <- addLog(mod, paste0("Model created from a ", src, " via fromlavaan()."))
  mod@log[[1]] <- NULL # Remove the auto "Model created" entry.

  if (run) mod <- runmodel(mod)
  mod
}
