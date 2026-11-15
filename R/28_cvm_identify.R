# Rows of the parameter table that can never be estimated in a cvm model:
#   zero:      structural zeros (factor loadings of composite indicators or on
#              composites, weights outside the composite of an indicator, and
#              elements of the "_c" matrices connecting different composites);
#   intercept: the intercepts of composites (computed from the indicator means);
#   variance:  the (residual) variances or scalings of composites (computed).
cvm_structural_rows <- function(x, pars = x@parameters){
  comp <- as.vector(as.matrix(x@extramatrices$compvec)) == 1
  block <- as.vector(as.matrix(x@extramatrices$blockvec))
  compIdx <- which(comp)
  varMat <- switch(x@types$latent, "cov" = "sigma_zeta", "ggm" = "delta_zeta")
  crossMats <- c("sigma_c", "kappa_c", "omega_c", "lowertri_c", "rho_c", "delta_c", "SD_c")

  isLambda <- pars$matrix == "lambda"
  isWeights <- pars$matrix == "weights"
  isCross <- pars$matrix %in% crossMats

  zero <- rep(FALSE, nrow(pars))
  zero[isLambda] <- block[pars$row[isLambda]] != 0 | pars$col[isLambda] %in% compIdx
  zero[isWeights] <- block[pars$row[isWeights]] != pars$col[isWeights]
  zero[isCross] <- block[pars$var1_id[isCross]] != block[pars$var2_id[isCross]]

  list(
    zero = zero,
    intercept = pars$matrix == "nu_eta" & pars$row %in% compIdx,
    variance = pars$matrix == varMat & pars$row == pars$col & pars$row %in% compIdx
  )
}

# Identification of composite variable models (cvm).
#
# Next to the usual identification of the common factors (as in lvm), this
# function marks all entries that may never be estimated in a cvm model as
# fixed and 'identified' (so that modification indices and stepup() skip them):
#   - the (residual) variances / scalings and intercepts of composites, which
#     are computed from the weights and the composite indicator covariances;
#   - factor loadings of composite indicators and on composites;
#   - weights outside of the composite an indicator belongs to;
#   - covariances between indicators of different composites in the "_c"
#     matrices (these are implied by the relations between the composites).
identify_cvm <- function(x, identify = TRUE){
  pars <- x@parameters
  type <- x@identification
  latent <- x@types$latent
  nGroup <- nrow(x@sample@groups)

  # Which constructs are composites, and to which composite do indicators belong?
  comp <- as.vector(as.matrix(x@extramatrices$compvec)) == 1
  block <- as.vector(as.matrix(x@extramatrices$blockvec))
  compIdx <- which(comp)
  facIdx <- which(!comp)

  # Helper to fix parameters:
  fixrows <- function(pars, rows, value = 0){
    if (is.logical(rows)) rows <- which(rows)
    if (length(rows) == 0) return(pars)
    if (!is.null(value)) pars$est[rows] <- value
    pars$par[rows] <- 0
    pars$fixed[rows] <- TRUE
    pars$identified[rows] <- TRUE
    clearpars(pars, rows)
  }

  # Always set beta diagonal to zero:
  pars <- fixrows(pars, pars$matrix == "beta" & pars$row == pars$col)

  ### Structural zeros and computed entries ###
  struc <- cvm_structural_rows(x, pars)
  pars <- fixrows(pars, struc$zero)
  pars <- fixrows(pars, struc$intercept, value = NULL)
  pars <- fixrows(pars, struc$variance, value = NULL)
  # The computed entries keep the value from the last evaluation of the model (if
  # any); they are recomputed whenever the model is evaluated:
  novalue <- which(struc$variance & (!is.finite(pars$est) | pars$est <= 0))
  pars$est[novalue] <- 1

  ### Single indicator factors: residual (co)variances fixed to zero ###
  lamPars <- pars[pars$matrix == "lambda" & (!pars$fixed | pars$est != 0), , drop = FALSE]
  for (f in facIdx){
    inds <- unique(lamPars$var1_id[lamPars$col == f])
    if (length(inds) == 1){
      cons <- pars$var1_id == inds &
        pars$matrix %in% c("sigma_epsilon","omega_epsilon","delta_epsilon","lowertri_epsilon","kappa_epsilon","rho_epsilon","SD_epsilon")
      pars <- fixrows(pars, cons)
    }
  }

  # The remainder (scaling of composites and factors) is only done when the
  # model should be identified automatically:
  if (identify){
    ### Scaling of the composites: fix the first free weight to one ###
    for (g in seq_len(nGroup)){
      for (k in compIdx){
        inComp <- pars$matrix == "weights" & pars$col == k & pars$group_id == g
        # If a fixed nonzero weight already scales the composite, leave it as is:
        if (any(pars$fixed[inComp] & pars$est[inComp] != 0, na.rm = TRUE)) next
        first <- which(inComp & !pars$fixed)[1]
        if (is.na(first)) next
        pars <- fixrows(pars, first, value = 1)
      }
    }
  
    ### Identification of the common factors (as in lvm) ###
    facVarMat <- switch(latent, "cov" = "sigma_zeta", "ggm" = "delta_zeta")
    isFacMean <- pars$matrix == "nu_eta" & pars$row %in% facIdx
    isFacVar <- pars$matrix == facVarMat & pars$row == pars$col & pars$row %in% facIdx
    isFacNu <- pars$matrix == "nu" & block[pars$row] == 0
    isFacLambda <- pars$matrix == "lambda" & pars$col %in% facIdx & block[pars$row] == 0
  
    if (length(facIdx) > 0){
      if (nGroup == 1){
        # Latent means: fixed to zero unless enough intercepts are fixed:
        nFixedNu <- sum(pars$fixed[isFacNu])
        if (nFixedNu <= length(facIdx)){
          pars <- fixrows(pars, isFacMean)
        }
  
        if (type == "variance"){
          pars <- fixrows(pars, isFacVar, value = 1)
        } else {
          for (f in facIdx){
            inFac <- isFacLambda & pars$col == f
            if (any(pars$fixed[inFac] & pars$est[inFac] != 0, na.rm = TRUE)) next
            first <- which(inFac & !pars$fixed)[1]
            if (is.na(first)) next
            pars <- fixrows(pars, first, value = 1)
          }
        }
      } else {
        # Number of equality constrained parameters per element:
        eqCount <- function(sel){
          if (!any(sel)) return(0)
          sub <- pars[sel, , drop = FALSE]
          key <- paste(sub$row, sub$col)
          sum(tapply(seq_len(nrow(sub)), key, function(i){
            !all(sub$fixed[i]) && length(unique(sub$par[i])) == 1
          }))
        }
  
        # Latent means: free in all but the first group if enough intercepts are equal:
        anyNuFixed <- any(pars$fixed[isFacNu] & !pars$identified[isFacNu])
        if (!anyNuFixed){
          if (eqCount(isFacNu) >= length(facIdx)){
            fixed <- isFacMean & pars$group_id == 1
            free <- which(isFacMean & pars$group_id > 1 & !(pars$fixed & !pars$identified))
          } else {
            fixed <- isFacMean
            free <- integer(0)
          }
          pars <- fixrows(pars, fixed)
          if (length(free) > 0){
            pars$par[free] <- max(pars$par) + seq_along(free)
            pars$fixed[free] <- FALSE
            pars$identified[free] <- FALSE
            pars <- clearpars(pars, free)
          }
        }
  
        if (type == "variance"){
          if (eqCount(isFacLambda) >= length(facIdx)){
            fixed <- isFacVar & pars$group_id == 1
            free <- which(isFacVar & pars$group_id > 1 & pars$identified)
          } else {
            fixed <- isFacVar
            free <- integer(0)
          }
          pars <- fixrows(pars, fixed, value = 1)
          if (length(free) > 0){
            pars$par[free] <- max(pars$par) + seq_along(free)
            pars$fixed[free] <- FALSE
            pars$identified[free] <- FALSE
            pars <- clearpars(pars, free)
          }
        } else {
          for (g in seq_len(nGroup)){
            for (f in facIdx){
              inFac <- isFacLambda & pars$col == f & pars$group_id == g
              if (any(pars$fixed[inFac] & pars$est[inFac] != 0, na.rm = TRUE)) next
              first <- which(inFac & !pars$fixed)[1]
              if (is.na(first)) next
              pars <- fixrows(pars, first, value = 1)
            }
          }
        }
      }
    }
  
  }

  # Fix labels:
  x@parameters <- parRelabel(pars)

  # Update which structural parameters may be freed (latent = "ggm"):
  x <- cvm_update_structure(x)

  # Return model:
  return(x)
}

# Which structural parameters of a cvm model with a latent network may be freed.
# This is part of the identification step: it is only called by identify_cvm().
#
# With latent = "ggm" the scaling of a composite (delta_zeta) is computed from
# its variance. This is only possible in closed form when a composite either
# (a) has no incoming effects in beta, in which case it can be a node in the
# latent network (omega_zeta), or (b) has incoming effects but no edges in the
# latent network. A composite can therefore not have both; the implied model
# (implied_cvm, R and C++) stops with an error for such a model. This function
# marks the parameters that would lead to such a model as 'identified', so that
# modification indices and model search skip them:
#  - for a composite with incoming effects: its edges in the latent network;
#  - for a composite without incoming effects (a node of the latent network,
#    whether or not it currently has edges): its incoming effects.
# The marks only depend on which composites have incoming effects. Model search
# does not free marked parameters and can therefore not change that, so the
# marks stay valid during stepup() without updating them; prune() and the other
# functions that remove parameters identify the model again afterwards. If
# effects are added or removed by hand (freepar(), fixpar()), the marks are
# updated the next time the model is identified (identify()).
cvm_update_structure <- function(x){
  pars <- x@parameters
  comp <- which(as.vector(as.matrix(x@extramatrices$compvec)) == 1)
  if (length(comp) == 0) return(x)

  if (!identical(x@types$latent, "ggm")){
    # No restriction when the constructs are modeled with covariances: remove
    # marks on effects that were set while they formed a network (the type of
    # a model can be changed with transmod()):
    inactiveBeta <- pars$matrix == "beta" & pars$row != pars$col & pars$row %in% comp &
      pars$fixed & !is.na(pars$est) & pars$est == 0
    pars$identified[inactiveBeta] <- FALSE
    x@parameters <- pars
    return(x)
  }

  nLat <- length(as.vector(as.matrix(x@extramatrices$compvec)))
  latents <- unique(pars$var1[pars$matrix == "nu_eta"])

  # Parameters that are in the model (free, or fixed to a nonzero value):
  active <- !pars$fixed | (!is.na(pars$est) & pars$est != 0)

  for (g in unique(pars$group_id)){
    isBeta <- pars$matrix == "beta" & pars$group_id == g & pars$row != pars$col
    isOmega <- pars$matrix == "omega_zeta" & pars$group_id == g

    lockBeta <- rep(FALSE, nLat)
    lockOmega <- rep(FALSE, nLat)
    for (k in comp){
      hasBeta <- any(active[isBeta & pars$row == k])
      hasOmega <- any(active[isOmega & (pars$row == k | pars$col == k)])
      if (hasBeta && hasOmega){
        stop(paste0("Composite '", latents[k], "' has both incoming effects in 'beta' and edges in the latent network 'omega_zeta'. ",
                    "With latent = 'ggm', a composite can only be part of the latent network if it has no incoming effects ",
                    "(its scaling is computed from the composite variance). Remove either the effects or the edges, or use latent = 'cov'."))
      }
      lockOmega[k] <- hasBeta
      lockBeta[k] <- !hasBeta
    }

    # Effects on composites that are nodes of the latent network:
    inactiveBeta <- which(isBeta & !active & pars$row %in% comp)
    pars$identified[inactiveBeta] <- lockBeta[pars$row[inactiveBeta]]

    # Edges of composites that have incoming effects:
    inactiveOmega <- which(isOmega & !active & (pars$row %in% comp | pars$col %in% comp))
    pars$identified[inactiveOmega] <- lockOmega[pars$row[inactiveOmega]] | lockOmega[pars$col[inactiveOmega]]
  }

  x@parameters <- pars
  x
}

# Support for transmod() in the cvm family: checks the requested types, and
# returns which rows of the parameter table are structural elements, i.e.,
# elements that are never estimated and are therefore not part of the structure
# that transmod() transforms:
#  - elements between indicators of different composites in the "_c" matrices
#    (these matrices are block diagonal with one block per composite, and every
#    parameterization of a block diagonal covariance matrix is block diagonal);
#  - loadings of composite indicators and weights outside of a composite.
# The computed (residual) variances / scalings of composites are part of the
# structure: their values are used to transform the construct covariances.
cvm_transmod_setup <- function(x, dots){
  if ("latent" %in% names(dots) && !dots[["latent"]] %in% c("cov", "ggm")){
    stop("In the cvm family, 'latent' can only be 'cov' or 'ggm'.")
  }
  for (type in intersect(names(dots), c("residual", "composite"))){
    if (!dots[[type]] %in% c("cov", "chol", "prec", "ggm", "cor")){
      stop(paste0("'", type, "' must be one of 'cov', 'chol', 'prec', 'ggm' or 'cor'."))
    }
  }
  cvm_structural_rows(x)$zero
}

# Default matrices for model search (stepup, prune, modelsearch, ...) in cvm
# models: the network structures that are in the model.
cvm_default_matrices <- function(x){
  matrices <- character(0)
  if (identical(x@types$latent, "ggm")){
    matrices <- c(matrices, "omega_zeta")
  }
  if (identical(x@types$residual, "ggm")){
    matrices <- c(matrices, "omega_epsilon")
  } else if (identical(x@types$residual, "prec")){
    matrices <- c(matrices, "kappa_epsilon")
  }
  if (identical(x@types$composite, "ggm")){
    matrices <- c(matrices, "omega_c")
  } else if (identical(x@types$composite, "prec")){
    matrices <- c(matrices, "kappa_c")
  }
  matrices[matrices %in% x@matrices$name]
}

# Write the computed quantities of the composites (their (residual) variance or
# scaling and their intercept) to the parameter table, so that these are shown
# next to the estimated parameters. These entries remain fixed: they are
# recomputed from the weights and the composite indicator covariances every
# time the model is evaluated, regardless of the value in the parameter table.
cvm_fill_computed <- function(x){
  comp <- which(as.vector(as.matrix(x@extramatrices$compvec)) == 1)
  if (length(comp) == 0 || length(x@modelmatrices) == 0) return(x)

  pars <- x@parameters
  varMat <- switch(x@types$latent, "cov" = "sigma_zeta", "ggm" = "delta_zeta")

  for (g in seq_along(x@modelmatrices)){
    mats <- x@modelmatrices[[g]]
    if (!is.null(mats[[varMat]])){
      rows <- which(pars$group_id == g & pars$matrix == varMat & pars$row == pars$col & pars$row %in% comp)
      pars$est[rows] <- diag(as.matrix(mats[[varMat]]))[pars$row[rows]]
    }
    if (!is.null(mats[["nu_eta"]])){
      rows <- which(pars$group_id == g & pars$matrix == "nu_eta" & pars$row %in% comp)
      pars$est[rows] <- as.vector(as.matrix(mats[["nu_eta"]]))[pars$row[rows]]
    }
  }

  # Warn when a computed residual variance of a composite is not positive:
  latents <- unique(pars$var1[pars$matrix == "nu_eta"])
  negative <- character(0)
  for (g in seq_along(x@modelmatrices)){
    sz <- x@modelmatrices[[g]][["sigma_zeta"]]
    if (!is.null(sz)){
      neg <- comp[!(diag(as.matrix(sz))[comp] > 0)]
      if (length(neg) > 0){
        negative <- c(negative, if (length(latents) >= max(neg)) latents[neg] else as.character(neg))
      }
    }
  }
  if (length(negative) > 0){
    warning(paste0("The computed residual variance of the following composite(s) is not positive: ", paste(unique(negative), collapse = ", "),
                   ". The variance explained by the predictors of a composite then exceeds the variance of the composite itself ",
                   "(an improper solution); interpret the results with care."), call. = FALSE)
  }

  x@parameters <- pars
  x
}

# Composites with estimated weights that are not related to any other
# construct: their weights are not identified.
cvm_unconnected_composites <- function(x){
  pars <- x@parameters
  comp <- which(as.vector(as.matrix(x@extramatrices$compvec)) == 1)
  latents <- unique(pars$var1[pars$matrix == "nu_eta"])
  if (length(latents) == 0) latents <- paste0("Eta_", seq_along(as.vector(as.matrix(x@extramatrices$compvec))))
  active <- !pars$fixed | (!is.na(pars$est) & pars$est != 0)
  structural <- pars$matrix %in% c("beta", "sigma_zeta", "omega_zeta") & pars$row != pars$col
  res <- character(0)
  for (g in unique(pars$group_id)){
    for (k in comp){
      freeWeights <- any(!pars$fixed[pars$matrix == "weights" & pars$col == k & pars$group_id == g])
      connected <- any(active[structural & pars$group_id == g & (pars$row == k | pars$col == k)])
      if (freeWeights && !connected) res <- c(res, latents[k])
    }
  }
  unique(res)
}

# Make sure the starting values of a cvm model are admissible: all computed
# (residual) variances of composites positive and all implied covariance
# matrices positive definite. Starting values are formed per group and per
# matrix, which, in particular in combination with equality constraints across
# groups, can lead to a composite of which the starting values of the effects
# explain more than its variance. In that case the starting values of the
# structural part (effects, and covariances / edges between constructs) are
# shrunk towards zero until the starting point is admissible (with all of these
# at zero the model is always admissible).
cvm_admissible_start <- function(x, maxtry = 25){
  # Make the estimates of equality constrained parameters consistent:
  x <- updateModel(parVector(x), x, updateMatrices = FALSE)

  admissible <- function(x){
    imp <- tryCatch(implied_cvm(x, all = FALSE), error = function(e) NULL)
    if (is.null(imp)) return(FALSE)
    all(vapply(imp, function(g){
      if (!isTRUE(g$proper_composites)) return(FALSE)
      sig <- as.matrix(g$sigma)
      if (any(!is.finite(sig))) return(FALSE)
      ev <- eigen(sig, symmetric = TRUE, only.values = TRUE)$values
      min(ev) > sqrt(.Machine$double.eps) * max(abs(ev))
    }, logical(1)))
  }

  pars <- x@parameters
  structural <- which(!pars$fixed & (pars$matrix == "beta" |
                                       (pars$matrix %in% c("sigma_zeta", "omega_zeta") & pars$row != pars$col)))
  if (length(structural) == 0) return(x)

  try <- 0
  while (try < maxtry && !admissible(x)){
    x@parameters$est[structural] <- 0.5 * x@parameters$est[structural]
    try <- try + 1
  }

  x
}

# Starting values for the intercepts of the indicators of common factors.
#
# The mean of a composite is the weighted sum of the means of its indicators,
# which is in general not zero. A common factor that is predicted by composites
# (or by factors that are) then has a nonzero implied mean as well, so that the
# intercepts of its indicators (nu) are the means of these indicators MINUS the
# loadings times that mean. The generic starting values put all intercepts at
# the sample means. For such a factor that can be far off (several units in
# the tests), so that estimation starts with a large misfit in the means; with
# the DWLS estimator this made the optimizer stop without convergence. This
# function starts the free intercepts of the factor indicators at
#   sample mean - lambda %*% (implied mean of the constructs)
# evaluated at the other starting values. The intercepts ("means") of composite
# indicators are not changed: these are the sample means themselves.
cvm_start_intercepts <- function(x){
  if (!isTRUE(x@meanstructure)) return(x)
  imp <- tryCatch(implied_cvm(x, all = TRUE), error = function(e) NULL)
  if (is.null(imp)) return(x)

  pars <- x@parameters
  block <- as.vector(as.matrix(x@extramatrices$blockvec))
  changed <- FALSE
  for (g in seq_along(imp)){
    if (is.null(imp[[g]]$mu_eta) || length(x@sample@means) < g) next
    mu_eta <- as.vector(as.matrix(imp[[g]]$mu_eta))
    means <- as.vector(as.matrix(x@sample@means[[g]]))
    if (any(!is.finite(mu_eta)) || any(!is.finite(means)) || length(means) != length(block)) next

    # Loadings of the common factors at the starting values:
    lam <- pars[pars$matrix == "lambda" & pars$group_id == g, , drop = FALSE]
    Lf <- matrix(0, length(means), length(mu_eta))
    Lf[cbind(lam$row, lam$col)] <- lam$est

    rows <- which(pars$matrix == "nu" & pars$group_id == g & !pars$fixed & block[pars$row] == 0)
    if (length(rows) == 0) next
    start <- (means - as.vector(Lf %*% mu_eta))[pars$row[rows]]
    ok <- is.finite(start)
    pars$est[rows[ok]] <- start[ok]
    changed <- TRUE
  }
  if (!changed) return(x)
  x@parameters <- pars

  # Make the estimates of equality constrained parameters consistent:
  updateModel(parVector(x), x, updateMatrices = FALSE)
}

# Expected Hessian of the fit function (twice the unit Fisher information).
#
# NOT USED by the package: kept for later investigation, see the note at the
# nlminb settings in runmodel(). Nothing in this function is specific to cvm
# models. It can be supplied to nlminb for any model, which then performs
# Fisher scoring within a trust region instead of quasi-Newton iterations:
#   setoptimizer(x, "nlminb", optim.args = list(hess = psychonetrics:::cvm_expected_hessian))
# The expected Hessian is not the Hessian of the fit function: the two only
# coincide when the model reproduces the sample moments, and are close when
# the model fits well and the parameters are near the solution.
cvm_expected_hessian <- function(x, model){
  model <- updateModel(x, model, updateMatrices = FALSE)
  H <- tryCatch({
    if (model@cpp){
      2 * as.matrix(psychonetrics_FisherInformation_cpp(model, FALSE, FALSE))
    } else {
      2 * as.matrix(psychonetrics_FisherInformation(model))
    }
  }, error = function(e) NULL)

  # Fall back to a (scaled) identity matrix if the information could not be computed:
  if (is.null(H) || any(!is.finite(H))){
    H <- diag(length(x))
  }
  H
}

# Warnings on the solution of a cvm model, given after estimation:
#  - the optimizer did not report convergence;
#  - the indicator that scales a composite (weight fixed to 1) has a small
#    weight relative to the other indicators: a standardized weight below 1/25
#    of the largest of the composite, the ratio at which the starting weights
#    of a composite are replaced by equal weights in cvm(). The other weights
#    are then large and poorly determined in that scaling. In an example in
#    which the first indicator was forced to scale the composite
#    (bad_marker_crit = 0; standardized weight 1.4% of the largest at the
#    solution) the optimizer stopped 6.4 chi-square points before the optimum;
#  - the indicator that scales a composite correlates negatively with that
#    composite: the composite is reversed relative to its scaling indicator
#    (the weight of that indicator has another sign than its loading), which
#    goes together with large weights and made the optimizer stop before the
#    optimum in examples;
#  - a composite with estimated weights is not related to any other construct
#    (its weights are then not identified);
#  - an element that is not a parameter of a cvm model was freed by hand (the
#    identification step fixes these elements, see identify_cvm()).
cvm_solution_warnings <- function(x){
  msgs <- character(0)

  # Convergence of the optimizer:
  conv <- x@optim$convergence
  if (!is.null(conv) && length(conv) == 1 && !is.na(conv) && is.numeric(conv) && conv != 0){
    msg <- x@optim$message
    msgs <- c(msgs, paste0("the optimizer did not report convergence",
                           if (!is.null(msg) && length(msg) == 1) paste0(" ('", msg, "')") else "",
                           " (consider other starting values or another optimizer, see setoptimizer())"))
  }

  # Small weight of the indicator that scales a composite, or a composite that
  # is reversed relative to that indicator:
  pars <- x@parameters
  comp <- which(as.vector(as.matrix(x@extramatrices$compvec)) == 1)
  block <- as.vector(as.matrix(x@extramatrices$blockvec))
  compInd <- which(block != 0)
  latents <- unique(pars$var1[pars$matrix == "nu_eta"])
  # Composites with estimated weights that are not related to any other
  # construct (their weights are arbitrary, so the scaling is not checked):
  unconnected <- cvm_unconnected_composites(x)
  bad <- character(0)
  reversed <- character(0)
  for (g in seq_along(x@modelmatrices)){
    W <- x@modelmatrices[[g]][["weights"]]
    Tc <- x@modelmatrices[[g]][["sigma_c"]]
    if (is.null(W) || is.null(Tc)) next
    W <- as.matrix(W)
    Tc <- as.matrix(Tc)
    sds <- numeric(length(block))
    sds[compInd] <- sqrt(pmax(diag(Tc), 0))
    for (k in comp){
      idx <- which(block == k)
      if (length(idx) < 2) next
      name <- if (length(latents) >= k) latents[k] else as.character(k)
      if (name %in% unconnected) next
      isMarker <- pars$matrix == "weights" & pars$col == k & pars$group_id == g & pars$fixed & pars$est != 0
      anyFree <- any(pars$matrix == "weights" & pars$col == k & pars$group_id == g & !pars$fixed)
      if (sum(isMarker) != 1 || !anyFree) next
      rowMarker <- pars$row[isMarker]
      s <- abs(W[idx, k]) * sds[idx]
      sMarker <- abs(W[rowMarker, k]) * sds[rowMarker]
      if (all(is.finite(s)) && max(s) > 0 && sMarker < max(s) / 25){
        bad <- c(bad, name)
      }
      # Covariance of the scaling indicator with the composite:
      idxc <- match(idx, compInd)
      covMarker <- sum(Tc[match(rowMarker, compInd), idxc] * W[idx, k])
      if (is.finite(covMarker) && sign(W[rowMarker, k]) * covMarker < 0){
        reversed <- c(reversed, name)
      }
    }
  }
  if (length(bad) > 0){
    msgs <- c(msgs, paste0("the indicator that scales composite(s) ", paste(unique(bad), collapse = ", "),
                           " (weight fixed to 1) has a small weight relative to the other indicators (a standardized weight below 1/25 of the largest), so that the other weights are large and poorly determined in this scaling and the optimizer can stop before the optimum",
                           " (consider scaling the composite with another indicator: free the fixed weight with freepar() and fix the weight of another indicator to 1 with fixpar())"))
  }
  if (length(reversed) > 0){
    msgs <- c(msgs, paste0("the indicator that scales composite(s) ", paste(unique(reversed), collapse = ", "),
                           " (weight fixed to 1) correlates negatively with that composite, so that the composite is reversed relative to its scaling indicator;",
                           " this does not change the fit of the model at the optimum, but the other weights tend to be large and the optimizer can stop before the optimum",
                           " (consider scaling the composite with another indicator: use the default 'bad_marker_crit = 1' in cvm(), or free the fixed weight with freepar() and fix the weight of another indicator to 1 with fixpar())"))
  }

  # Composites with estimated weights that are not related to any other construct:
  if (length(unconnected) > 0){
    msgs <- c(msgs, paste0("the following composite(s) have estimated weights but are not related to any other construct (no effects in 'beta' and no covariances / edges): ",
                           paste(unconnected, collapse = ", "), " (the weights of an isolated composite are not identified)"))
  }

  # Elements that were freed by hand but are not parameters of a cvm model:
  struc <- cvm_structural_rows(x)
  notpar <- which((struc$zero | struc$intercept | struc$variance) & !pars$fixed)
  if (length(notpar) > 0){
    msgs <- c(msgs, paste0("the following element(s) were freed but cannot be estimated in a composite variable model, because they are computed from other parameters ",
                           "(variances / scalings and intercepts of composites) or are structural zeros (loadings of composite indicators, weights outside ",
                           "the composite of an indicator, and elements connecting indicators of different composites in the '_c' matrices): ",
                           paste(unique(paste0(pars$matrix[notpar], "[", pars$var1[notpar], ", ",
                                               ifelse(is.na(pars$var2[notpar]), "1", pars$var2[notpar]), "]")), collapse = "; "),
                           " (fix these again with fixpar())"))
  }

  if (length(msgs) > 0){
    warning(paste0("The estimates of this composite variable model may not be reliable: ", paste(msgs, collapse = "; in addition, "), "."), call. = FALSE)
  }
  invisible(NULL)
}
