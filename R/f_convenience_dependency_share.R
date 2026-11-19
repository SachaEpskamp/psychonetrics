# The share of the dependence between the composites that is carried by a
# composite model:
#
#   share = 1 - (chisq - df) / (chisq_0 - df_0)
#
# in which chisq_0 and df_0 belong to a baseline model in which the indicators
# of different composites are independent, while all means, variances and
# covariances within composites are free. This is the relative noncentrality
# index (the CFI without truncation) with this baseline model in place of the
# default baseline model.

# Checks the model and returns the set of every indicator. The sets are the
# composites, in which all composites with a single indicator form one set: a
# model with a saturated structure between the composites reproduces the
# relations among these indicators exactly.
dependency_share_sets <- function(x){
  # 1. A cvm model:
  if (!is(x, "psychonetrics")){
    stop("'x' must be a psychonetrics model made with cvm() or cnm().")
  }
  if (!identical(x@model, "cvm")){
    stop("The share of the dependency is only available for composite models: 'x' is a model of the '", x@model,
         "' family. Use a model made with cvm() or cnm().")
  }

  # 2. Computed:
  if (!isTRUE(x@computed)){
    stop("Model has not yet been computed. Run the model first with runmodel().")
  }
  if (is.null(x@fitmeasures$chisq) || is.null(x@fitmeasures$df)){
    stop("No fit measures computed yet. Run the model with runmodel(), using the default 'addfit = TRUE'.")
  }

  # 3. A single group:
  if (nrow(x@sample@groups) != 1){
    stop("Models with multiple groups are not supported yet: the share of the dependency is only available for a single group.")
  }

  # 4. ML or FIML:
  if (!x@estimator %in% c("ML", "FIML")){
    stop("The share of the dependency is only available for the 'ML' and 'FIML' estimators, and this model uses the '",
         x@estimator, "' estimator. Estimate the model with estimator = \"ML\" or estimator = \"FIML\".")
  }

  pars <- x@parameters
  vars <- x@sample@variables$label
  nVar <- length(vars)
  isSet <- !pars$fixed | pars$est != 0 # Elements that are in the model

  # Pattern of the weights (variables x constructs):
  wt <- pars[pars$matrix == "weights", , drop = FALSE]
  nLatent <- max(c(pars$col[pars$matrix %in% c("weights", "lambda")], 0))
  if (nrow(wt) == 0 || nLatent == 0){
    stop("No composite weights found in the model.")
  }
  latents <- rep(NA_character_, nLatent)
  latents[wt$col] <- wt$var2
  W <- matrix(FALSE, nVar, nLatent)
  W[cbind(wt$row, wt$col)] <- isSet[pars$matrix == "weights"]

  # 5. All constructs are composites:
  lam <- pars[pars$matrix == "lambda" & isSet, , drop = FALSE]
  if (nrow(lam) > 0 || any(colSums(W) == 0)){
    notComp <- unique(c(lam$var2, latents[colSums(W) == 0]))
    stop("All constructs must be composites, and this model contains ", length(notComp), " common factor(s)",
         if (!all(is.na(notComp))) paste0(" (", paste(notComp[!is.na(notComp)], collapse = ", "), ")"),
         ". The share of the dependency requires a model with only 'weights' and no 'lambda'.")
  }

  # 6. Every indicator belongs to exactly one composite:
  nComp <- rowSums(W)
  if (any(nComp != 1)){
    stop("Every indicator must belong to exactly one composite, as given by its nonzero weights. ",
         if (any(nComp == 0)) paste0("Without a nonzero weight: ", paste(vars[nComp == 0], collapse = ", "), ". "),
         if (any(nComp > 1)) paste0("In more than one composite: ", paste(vars[nComp > 1], collapse = ", "), ". "),
         "Use a model in which every variable has a (free or fixed) nonzero weight on one composite.")
  }
  composite <- as.vector(W %*% seq_len(nLatent))

  # 7. No regressions between the composites:
  nBeta <- sum(pars$matrix == "beta" & isSet)
  if (nBeta > 0){
    stop("The model contains ", nBeta, " effect(s) between the composites in 'beta'. The share of the dependency requires ",
         "a model without 'beta', in which all composites are related through ",
         if (identical(x@types$latent, "ggm")) "edges ('omega_zeta')." else "covariances ('sigma_zeta').")
  }

  # Parameters that are equal to another parameter are not free of each other:
  freePar <- pars$par[!pars$fixed & pars$par != 0]
  constrained <- pars$par %in% freePar[duplicated(freePar)]

  # 8. Saturated structure between the composites:
  ggm <- identical(x@types$latent, "ggm")
  zetaMat <- if (ggm) "omega_zeta" else "sigma_zeta"
  what <- if (ggm) "edges" else "covariances"
  between <- pars$matrix == zetaMat & pars$row != pars$col
  nBetween <- nLatent * (nLatent - 1) / 2
  nFree <- sum(between & !pars$fixed)
  if (nFree != nBetween){
    stop("The structure between the composites is not saturated: ", nBetween - nFree, " of the ", nBetween, " ", what,
         " between the composites ('", zetaMat, "') are fixed. The share of the dependency requires the model with all ", what,
         " between the composites, e.g., not a model after prune() or a model with fixed ", what, ".")
  }
  if (any(between & constrained)){
    stop("The structure between the composites is not saturated: ", sum(between & constrained), " of the ", what,
         " between the composites ('", zetaMat, "') are constrained to be equal to another parameter. ",
         "The share of the dependency requires the model in which all these ", what, " are free.")
  }

  # 9. Saturated structure within every composite:
  type <- x@types$composite
  offdiagMat <- switch(type, cov = "sigma_c", chol = "lowertri_c", prec = "kappa_c", ggm = "omega_c", cor = "rho_c",
                       stop("Unknown type of the composite indicators: '", type, "'."))
  diagMat <- switch(type, ggm = "delta_c", cor = "SD_c", offdiagMat)
  within <- pars$matrix == offdiagMat & pars$row != pars$col
  within[within] <- composite[pars$row[within]] == composite[pars$col[within]]
  within <- within | (pars$matrix == diagMat & pars$row == pars$col)
  nWithin <- sum(vapply(split(seq_len(nVar), composite), function(i) length(i) * (length(i) + 1) / 2, numeric(1)))
  nFree <- sum(within & !pars$fixed)
  if (nFree != nWithin){
    stop("The structure within the composites is not saturated: ", nWithin - nFree, " of the ", nWithin,
         " elements within the composites ('", paste(unique(c(offdiagMat, diagMat)), collapse = "' and '"), "') are fixed. ",
         "The share of the dependency requires that all variances and covariances of the indicators within a composite are free; ",
         "only the weights may be fixed.")
  }
  if (any(within & constrained)){
    stop("The structure within the composites is not saturated: ", sum(within & constrained),
         " of the elements within the composites ('", paste(unique(c(offdiagMat, diagMat)), collapse = "' and '"),
         "') are constrained to be equal to another parameter. The share of the dependency requires that these are free.")
  }

  # The means are saturated as well:
  if (isTRUE(x@meanstructure)){
    nu <- pars$matrix == "nu"
    if (sum(nu & !pars$fixed) != nVar || any(nu & constrained)){
      stop("The share of the dependency requires that all intercepts ('nu') are free: the baseline model has free means.")
    }
  }

  # The sets: composites, with all single-indicator composites merged into one:
  nIndicators <- tabulate(composite, nLatent)
  set <- ifelse(nIndicators[composite] == 1, 0, composite)
  if (length(unique(set)) < 2){
    stop("The share of the dependency is not defined for this model: ",
         if (all(set == 0)) "all composites have a single indicator, so the model is saturated. "
         else "there is only one composite. ",
         "At least two composites are needed, of which at least one has more than one indicator.")
  }
  names(set) <- vars
  set
}

# The baseline model: indicators in different sets are independent, and all
# means, variances and covariances within sets are free. It uses the sample
# statistics of x, and with that the same observations and missing-data
# handling.
dependency_share_baseline <- function(x, set){
  vars <- x@sample@variables$label
  pattern <- 1 * outer(unname(set), unname(set), "==")
  dimnames(pattern) <- list(vars, vars)
  baseline <- varcov(type = "cov", sigma = pattern, vars = vars, estimator = x@estimator,
                     sampleStats = x@sample, meanstructure = isTRUE(x@meanstructure), verbose = FALSE)
  baseline <- runmodel(baseline, addMIs = FALSE, addSEs = FALSE, addInformation = FALSE, verbose = FALSE)

  # The degrees of freedom are the covariances between indicators of different sets:
  df0 <- sum(pattern[lower.tri(pattern)] == 0)
  if (!isTRUE(all.equal(baseline@fitmeasures$df, df0))){
    stop("The baseline model has ", baseline@fitmeasures$df, " degrees of freedom, where ", df0,
         " were expected. Please report this to https://github.com/SachaEpskamp/psychonetrics/issues")
  }

  # Both chi-squares must be relative to the same saturated model:
  satModel <- x@fitmeasures$unrestricted.logl
  satBaseline <- baseline@fitmeasures$unrestricted.logl
  if (!is.null(satModel) && !is.null(satBaseline) && is.finite(satModel) && is.finite(satBaseline) &&
      abs(satModel - satBaseline) > 1e-3 + 1e-8 * abs(satModel)){
    stop("The saturated log-likelihood of the baseline model (", format(satBaseline, digits = 10),
         ") differs from that of the model (", format(satModel, digits = 10),
         "), so the two chi-squares cannot be compared. Re-estimate the model with runmodel(), using the default 'saturated' argument.")
  }
  baseline
}

dependency_share <- function(x, verbose = TRUE){
  experimentalWarning("dependency_share()")

  # Check the model and form the sets of indicators:
  set <- dependency_share_sets(x)

  # Estimate the baseline model:
  baseline <- dependency_share_baseline(x, set)
  chisq <- x@fitmeasures$chisq
  df <- x@fitmeasures$df
  chisq0 <- baseline@fitmeasures$chisq
  df0 <- baseline@fitmeasures$df

  # 10. The composites must be related:
  if (!is.finite(chisq0) || chisq0 <= df0){
    stop("The share of the dependency is not defined: the chi-square of the baseline model (", round(chisq0, 2),
         ") is not larger than its degrees of freedom (", df0, "), so the composites are not clearly related.")
  }

  # The share:
  share <- 1 - (chisq - df) / (chisq0 - df0)

  if (isTRUE(verbose)){
    cat("The calculated share of the dependency equals: ", formatC(share, digits = 3, format = "f"), "\n", sep = "")
  }
  invisible(share)
}
