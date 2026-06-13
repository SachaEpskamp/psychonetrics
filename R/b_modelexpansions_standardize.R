# Standardized (std.all) solutions for the lvm and varcov model families.
#
# Implements the completely standardized solution (Rosseel, 2012, JSS,
# doi:10.18637/jss.v048.i02; lavaan::standardizedSolution(type = "std.all")):
# every parameter is rescaled by the implied standard deviations of the
# (observed or latent) variables indexed by its row/column, so that all
# variables have unit variance. Point estimates are obtained from the implied
# model matrices per group; standard errors by the delta method
# (numDeriv::jacobian of the standardization map times getVCOV(x)).
#
# Only meaningful for the Gaussian lvm and varcov families; for every other
# model framework the standardized column is left NA (see parameters()).

# Per-matrix standardization rule for std.all. For a parameter occupying
# element (r, c) of model matrix M, the standardized value is by default
#   std[r, c] = est[r, c] * sd_row[r]^a * sd_col[c]^b
# with sd_row / sd_col the implied SD of the variable indexed by the row / col
# (observed or latent), and (a, b) the exponents below. Matrices that are
# already scale-free (partial correlations omega; correlations rho) use a = b = 0
# (std = est). Returns NULL for a matrix that has no defined standardization.
#
#  rowtype / coltype: "obs" (observed), "lat" (latent) or NA (no variable on
#  that side, e.g. an intercept's column). diagonal_only = TRUE means only the
#  diagonal (r == c) carries the row exponent (the SD-scaling matrices delta/SD).
#
#  kind = "residual" marks the residual / disturbance covariance matrices
#  (sigma_epsilon, sigma_zeta), which lavaan std.all standardizes specially
#  (Rosseel 2012): the DIAGONAL (residual / disturbance variances) is divided by
#  the TOTAL implied variance of that variable (sd_row[i]^2), turning it into the
#  proportion of variance that is residual; the OFF-DIAGONAL (residual /
#  disturbance covariances) is turned into a residual correlation by dividing by
#  sqrt(M[i,i] * M[j,j]) (the matrix's own diagonal), NOT by the total variances.
std_matrix_rule <- function(mat){
  switch(mat,
    # --- lvm matrices ---
    "lambda"        = list(a = -1, b =  1, rowtype = "obs", coltype = "lat"),
    "beta"          = list(a = -1, b =  1, rowtype = "lat", coltype = "lat"),
    "sigma_zeta"    = list(rowtype = "lat", coltype = "lat", kind = "residual"),
    "lowertri_zeta" = list(a = -1, b =  0, rowtype = "lat", coltype = "lat"),
    "kappa_zeta"    = list(a =  1, b =  1, rowtype = "lat", coltype = "lat"),
    "omega_zeta"    = list(a =  0, b =  0, rowtype = "lat", coltype = "lat"),
    "delta_zeta"    = list(a = -1, b =  0, rowtype = "lat", coltype = "lat", diagonal_only = TRUE),
    "sigma_epsilon" = list(rowtype = "obs", coltype = "obs", kind = "residual"),
    "lowertri_epsilon" = list(a = -1, b =  0, rowtype = "obs", coltype = "obs"),
    "kappa_epsilon" = list(a =  1, b =  1, rowtype = "obs", coltype = "obs"),
    "omega_epsilon" = list(a =  0, b =  0, rowtype = "obs", coltype = "obs"),
    "delta_epsilon" = list(a = -1, b =  0, rowtype = "obs", coltype = "obs", diagonal_only = TRUE),
    "nu"            = list(a = -1, b =  0, rowtype = "obs", coltype = NA),
    "nu_eta"        = list(a = -1, b =  0, rowtype = "lat", coltype = NA),
    # --- varcov matrices ---
    "sigma"         = list(a = -1, b = -1, rowtype = "obs", coltype = "obs"),
    "lowertri"      = list(a = -1, b =  0, rowtype = "obs", coltype = "obs"),
    "kappa"         = list(a =  1, b =  1, rowtype = "obs", coltype = "obs"),
    "omega"         = list(a =  0, b =  0, rowtype = "obs", coltype = "obs"),
    "delta"         = list(a = -1, b =  0, rowtype = "obs", coltype = "obs", diagonal_only = TRUE),
    "rho"           = list(a =  0, b =  0, rowtype = "obs", coltype = "obs"),
    "SD"            = list(a = -1, b =  0, rowtype = "obs", coltype = "obs", diagonal_only = TRUE),
    "mu"            = list(a = -1, b =  0, rowtype = "obs", coltype = NA),
    NULL)
}

# Per-group implied observed and latent standard deviations for the lvm /
# varcov families, from the implied model matrices in 'imp[[g]]'.
#   observed SDs: sqrt(diag(implied Sigma))
#   latent   SDs (lvm only): sqrt of the diagonal of the TOTAL (marginal) latent
#     covariance Cov(eta) = (I - B)^{-1} Sigma_zeta (I - B)^{-T}. This is the
#     marginal latent variance used to standardize loadings and regressions; for
#     a variance-identified model the total latent variance is 1 (so loadings are
#     scaled by sd_y alone), while for a regressed latent it exceeds the
#     disturbance variance Sigma_zeta[i, i].
std_implied_sds <- function(impg, family){
  sigma <- as.matrix(impg$sigma)
  sd_obs <- sqrt(diag(sigma))
  sd_lat <- NULL
  if (family == "lvm"){
    sigma_zeta <- as.matrix(impg$sigma_zeta)
    beta <- if (is.null(impg$beta)) matrix(0, nrow(sigma_zeta), nrow(sigma_zeta)) else as.matrix(impg$beta)
    nLat <- nrow(sigma_zeta)
    if (all(beta == 0)){
      eta_cov <- sigma_zeta
    } else {
      IBinv <- solve(diag(nLat) - beta)
      eta_cov <- IBinv %*% sigma_zeta %*% t(IBinv)
    }
    sd_lat <- sqrt(diag(eta_cov))
  }
  list(obs = sd_obs, lat = sd_lat)
}

# Standardize a single parameter-table value 'est' occupying element (row, col)
# of matrix 'mat' in group g, given the implied SDs (sds$obs, sds$lat) and, for
# the residual / disturbance covariance matrices, the implied diagonal of the
# matrix itself (mat_diag, the residual / disturbance variances). Returns 'est'
# unchanged for a matrix without a standardization rule.
std_one_value <- function(est, mat, row, col, sds, mat_diag = NULL){
  rule <- std_matrix_rule(mat)
  if (is.null(rule)) return(est)
  sd_row_vec <- if (identical(rule$rowtype, "lat")) sds$lat else sds$obs
  sd_col_vec <- if (identical(rule$coltype, "lat")) sds$lat else sds$obs

  # Residual / disturbance covariance matrices (sigma_epsilon, sigma_zeta):
  if (identical(rule$kind, "residual")){
    if (row == col){
      # Residual / disturbance variance -> proportion of TOTAL implied variance:
      return(est / sd_row_vec[row]^2)
    } else {
      # Residual / disturbance covariance -> residual correlation (within-matrix):
      if (is.null(mat_diag)) return(est)
      denom <- sqrt(mat_diag[row] * mat_diag[col])
      if (!is.finite(denom) || denom == 0) return(est)
      return(est / denom)
    }
  }

  if (isTRUE(rule$diagonal_only) && row != col) return(est)
  fac <- 1
  if (!is.null(rule$a) && rule$a != 0) fac <- fac * sd_row_vec[row]^rule$a
  if (!is.na(rule$coltype) && !is.null(rule$b) && rule$b != 0) fac <- fac * sd_col_vec[col]^rule$b
  est * fac
}

# Compute the std.all standardized value for every row of the parameter table
# of a fitted lvm / varcov model, from a free-parameter vector 'theta'. Used both
# for the point estimates (theta = parVector(x)) and as the map differentiated by
# the delta method. Returns a numeric vector with one entry per parameter-table
# row (NA for rows whose matrix has no standardization rule, e.g. unused).
std_values_per_row <- function(x, theta){
  partable <- x@parameters
  # Insert theta into the (first-row) free parameters and propagate to all rows
  # sharing the same par index (equality constraints):
  if (length(theta) > 0){
    free <- partable$par != 0
    partable$est[free] <- theta[partable$par[free]]
  }
  xt <- x
  xt@parameters <- partable
  imp <- impliedModel(xt, xt@types, all = FALSE)
  family <- x@model
  sds_per_group <- lapply(imp, std_implied_sds, family = family)

  out <- rep(NA_real_, nrow(partable))
  gids <- partable$group_id
  # Map group_id -> index in imp (imp is ordered by group order):
  group_order <- x@sample@groups$id
  for (gi in seq_along(imp)){
    gid <- group_order[gi]
    sds <- sds_per_group[[gi]]
    rows <- which(gids == gid)
    # Cache the diagonal of each residual / disturbance covariance matrix in this
    # group (needed to standardize their off-diagonal entries as correlations):
    diag_cache <- list()
    for (r in rows){
      mat <- partable$matrix[r]
      rule <- std_matrix_rule(mat)
      mat_diag <- NULL
      if (!is.null(rule) && identical(rule$kind, "residual")){
        if (is.null(diag_cache[[mat]]) && !is.null(imp[[gi]][[mat]])){
          diag_cache[[mat]] <- diag(as.matrix(imp[[gi]][[mat]]))
        }
        mat_diag <- diag_cache[[mat]]
      }
      out[r] <- std_one_value(partable$est[r], mat,
                              partable$row[r], partable$col[r], sds, mat_diag)
    }
  }
  out
}

# Populate the std (and, if SEs are available, se_std) columns of a fitted
# model's parameter table with the completely standardized (std.all) solution.
# Returns the model with x@parameters$std / x@parameters$se_std filled in.
# Only the lvm and varcov families are supported; the caller guards on that.
addStandardized <- function(x){
  if (!x@model %in% c("lvm", "varcov")){
    return(x)
  }
  theta <- parVector(x)

  # Point estimates: standardized value per parameter-table row.
  std_row <- tryCatch(std_values_per_row(x, theta), error = function(e) NULL)
  if (is.null(std_row)) return(x)
  x@parameters$std <- std_row

  # Ensure the se_std column exists (guarded for objects created before it was
  # added to the prototype):
  if (is.null(x@parameters[["se_std"]])){
    x@parameters$se_std <- rep(NA_real_, nrow(x@parameters))
  } else {
    x@parameters$se_std[] <- NA_real_
  }

  # Standard errors by the delta method, only if the model is computed and has
  # free parameters:
  if (!x@computed || length(theta) == 0 || all(x@parameters$par == 0)){
    return(x)
  }

  VCOV <- tryCatch(getVCOV(x), error = function(e) NULL)
  if (is.null(VCOV) || any(!is.finite(VCOV))) return(x)

  # The standardization map collapsed to one value per free parameter (the
  # first parameter-table row for each par index), so the Jacobian is
  # npar x npar and aligns with VCOV.
  fpr <- firstParRows(x)
  std_map <- function(par_vec){
    sv <- std_values_per_row(x, par_vec)
    sv[fpr]
  }
  J <- tryCatch(numDeriv::jacobian(std_map, theta), error = function(e) NULL)
  if (is.null(J) || any(!is.finite(J))) return(x)

  std_vcov <- J %*% VCOV %*% t(J)
  se_par <- sqrt(pmax(diag(std_vcov), 0))

  # Broadcast the per-par SE back to every parameter-table row (rows sharing a
  # par index get the same SE; fixed rows keep NA):
  free <- x@parameters$par != 0
  x@parameters$se_std[free] <- se_par[x@parameters$par[free]]

  x
}
