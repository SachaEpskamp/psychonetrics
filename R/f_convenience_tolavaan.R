# Convert a psychonetrics lvm (latent = "cov", residual = "cov") or cvm (in
# addition composite = "cov") into a fitted lavaan object or a lavaan parameter
# table.
#
# Only the lavaan public API (lavaan::lavaan) is used; no lavaan source code is
# ported. See man/tolavaan.Rd for the restriction to cov/cov models and the
# warm-start / storedata behavior.
tolavaan <- function(x, fit = TRUE, verbose = FALSE, ...){
  if (verbose) experimentalWarning("tolavaan()")

  ## ---- Guards (informative, in order) -------------------------------
  if (!is(x, "psychonetrics")){
    stop("'x' must be a psychonetrics object.")
  }
  if (!x@model %in% c("lvm", "cvm")){
    stop("Only lvm-family and cvm-family models can be converted to lavaan (got model = '", x@model, "').")
  }
  if (x@types$latent != "cov" || x@types$residual != "cov"){
    stop("Only latent = 'cov' and residual = 'cov' models can be converted; got latent = '",
         x@types$latent, "', residual = '", x@types$residual,
         "'. Network (ggm/prec) and Cholesky structures have no lavaan equivalent.")
  }
  isCVM <- x@model == "cvm"
  if (isCVM){
    if (x@types$composite != "cov"){
      stop("Only composite = 'cov' models can be converted; got composite = '", x@types$composite,
           "'. Network (ggm/prec), Cholesky and correlation structures of the composite indicators have no lavaan equivalent.")
    }
    if (utils::packageVersion("lavaan") < "0.7.0"){
      stop("Converting a composite variable model requires lavaan >= 0.7.")
    }
  }
  if (isTRUE(x@sample@corinput) || any(x@sample@variables$ordered)){
    stop("Ordinal / correlation-input models cannot be converted to lavaan.")
  }
  if (!x@estimator %in% c("ML", "FIML")){
    stop("Only estimator = 'ML' or 'FIML' can be converted to lavaan (got '", x@estimator, "').")
  }
  if (x@estimator == "FIML" && nrow(x@sample@rawdata) == 0){
    stop("FIML conversion requires the raw data; re-fit the lvm with storedata = TRUE.")
  }

  P    <- x@parameters
  vars <- x@sample@variables$label
  # Latent names: lambda columns (group 1) ordered by col:
  lam     <- P[P$matrix == "lambda" & P$group_id == 1, ]
  latents <- unique(lam$var2[order(lam$col)])
  nG      <- nrow(x@sample@groups)

  ## ---- Rows of the lavaan parameter table ---------------------------
  # Emit a row for EVERY managed parameter-table row (including fixed-at-zero
  # rows) so lavaan does not add any defaults of its own. The residual and
  # composite indicator blocks of a cvm are indexed within their own block, so
  # the variable ids (v1, v2) are used for these.
  maprow <- function(m, r, c, v1, v2){
    switch(m,
      lambda        = c(latents[c], "=~", vars[r]),
      weights       = c(latents[c], "<~", vars[r]),
      sigma_zeta    = c(latents[c], "~~", latents[r]),
      sigma_epsilon = c(vars[v2],   "~~", vars[v1]),
      sigma_c       = c(vars[v2],   "~~", vars[v1]),
      beta          = c(latents[r], "~",  latents[c]),
      nu            = c(vars[r],    "~1", ""),
      nu_eta        = c(latents[r], "~1", ""),
      NULL)
  }
  keep <- P$matrix %in% c("lambda", "weights", "sigma_zeta", "sigma_epsilon", "sigma_c", "beta", "nu", "nu_eta")
  P <- P[keep, , drop = FALSE]

  # Composite variable models: the structural zeros (loadings of composite
  # indicators and on composites, weights outside the composite of an
  # indicator, covariances between indicators of different composites) are not
  # part of the lavaan model, and the variances and intercepts of composites
  # are computed by lavaan as well (fixed, without a value):
  computed <- rep(FALSE, nrow(P))
  if (isCVM){
    struc <- cvm_structural_rows(x, P)
    P <- P[!struc$zero, , drop = FALSE]
    computed <- (struc$intercept | struc$variance)[!struc$zero]
  }
  lor <- t(mapply(maprow, P$matrix, P$row, P$col, P$var1_id, P$var2_id))

  PT <- data.frame(lhs = lor[, 1], op = lor[, 2], rhs = lor[, 3],
                   block = P$group_id, group = P$group_id,
                   free = 0L, ustart = NA_real_, label = "",
                   plabel = paste0(".p", seq_len(nrow(P)), "."),
                   stringsAsFactors = FALSE)
  isfree <- !P$fixed
  PT$free[isfree] <- seq_len(sum(isfree))
  PT$ustart[!isfree] <- P$est[!isfree]                 # fixed values
  PT$ustart[computed] <- NA_real_                      # computed by lavaan (composites)
  if (x@computed) PT$ustart[isfree] <- P$est[isfree]   # warm starts on free rows

  ## ---- Equality constraints: shared par index -> "==" rows ----------
  # lavaan does NOT enforce equality through duplicate free indices nor shared
  # labels; the only working encoding is per-row plabels plus extra "==" rows.
  eqrows <- NULL
  for (k in unique(P$par[P$par != 0])){
    w <- which(P$par == k)
    if (length(w) > 1){
      eqrows <- rbind(eqrows, data.frame(
        lhs = PT$plabel[w[1]], op = "==", rhs = PT$plabel[w[-1]],
        block = 0L, group = 0L, free = 0L, ustart = NA_real_,
        label = "", plabel = "", stringsAsFactors = FALSE))
    }
  }
  PT <- rbind(PT, eqrows)
  PT$id <- seq_len(nrow(PT))

  if (!fit) return(PT)

  ## ---- Data / sample statistics -------------------------------------
  lavargs <- list(model = PT, ...)
  glabels <- x@sample@groups$label
  if (nrow(x@sample@rawdata) > 0){                     # storedata = TRUE
    lavargs$data <- x@sample@rawdata
    if (nG > 1) lavargs$group <- x@sample@groupvar
    if (x@estimator == "FIML") lavargs$missing <- "ml"
  } else {
    covs <- lapply(x@sample@covs, as.matrix)
    for (g in seq_along(covs)) dimnames(covs[[g]]) <- list(vars, vars)
    names(covs) <- glabels
    lavargs$sample.cov  <- if (nG == 1) covs[[1]] else covs
    lavargs$sample.nobs <- x@sample@groups$nobs
    lavargs$sample.cov.rescale <- FALSE                # covs are already ML (n)
    lavargs$likelihood  <- "normal"
    if (x@meanstructure){
      mns <- lapply(x@sample@means, as.numeric)
      for (g in seq_along(mns)) names(mns[[g]]) <- vars
      names(mns) <- glabels
      lavargs$sample.mean <- if (nG == 1) mns[[1]] else mns
    }
  }
  do.call(lavaan::lavaan, lavargs)
}
