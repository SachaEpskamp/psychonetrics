# fixed.x: exogenous covariates whose means and mutual (co)variances are fixed to
# their sample values (and excluded from the free-parameter / statistic counts),
# matching lavaan's sem(..., fixed.x = TRUE) (Rosseel, 2012). Conditioning the
# model on x fixes the pure x-block at its saturated ML solution, which leaves
# every other ML estimate unchanged (so no estimator change is required); the
# x-block simply leaves the npar and the degrees-of-freedom statistic count.
#
# Two helpers below are shared by the lvm and varcov constructors:
#   apply_fixed_x_partable() : fix the x-block rows of a parameter table to the
#                              supplied sample moments and renumber the free pars.
#   fixed_x_nstat_drop()     : the number of sample statistics removed by fixing
#                              the x-block ( p_x means + p_x(p_x+1)/2 (co)variances
#                              per group ), used to reduce model@sample@nobs.

# Number of statistics removed from the df bookkeeping by fixing the x-block:
#   per group: p_x (x means, only when a mean structure is present) +
#              p_x(p_x+1)/2 (x variances + x-x covariances).
fixed_x_nstat_drop <- function(p_x, nGroup, meanstructure){
  if (p_x <= 0) return(0)
  vech_x <- p_x * (p_x + 1) / 2
  means_x <- if (isTRUE(meanstructure)) p_x else 0
  nGroup * (vech_x + means_x)
}

# Fix the x-block of a parameter table to the sample moments. 'x_idx' are the
# 1-based indices (into the variable ordering) of the exogenous variables.
# 'cov_matrix' is the name of the (symmetric) covariance matrix whose x-x block
# is fixed (e.g. "sigma" for varcov, "sigma_epsilon" for lvm), 'mean_matrix' the
# name of the mean/intercept matrix whose x rows are fixed ("mu" / "nu");
# mean_matrix may be NA when there is no mean structure. 'sample_covs' /
# 'sample_means' are per-group lists of the sample covariance matrices / mean
# vectors (in the variable ordering). Rows are fixed (fixed = TRUE, par = 0) and
# their est set to the corresponding sample moment; the free par indices are then
# renumbered to be contiguous. Returns the modified parameter table.
apply_fixed_x_partable <- function(partable, x_idx, cov_matrix, mean_matrix,
                                   sample_covs, sample_means, group_ids){
  if (length(x_idx) == 0) return(partable)

  for (gi in seq_along(group_ids)){
    gid <- group_ids[gi]
    S <- as.matrix(sample_covs[[gi]])
    grows <- partable$group_id == gid

    # Covariance x-x block: both row and col are x variables:
    csel <- grows & partable$matrix == cov_matrix &
      partable$row %in% x_idx & partable$col %in% x_idx
    if (any(csel)){
      idx <- which(csel)
      partable$est[idx] <- S[cbind(partable$row[idx], partable$col[idx])]
      partable$fixed[idx] <- TRUE
      partable$par[idx] <- 0
    }

    # Means / intercepts of the x variables:
    if (!is.na(mean_matrix) && !is.null(sample_means)){
      m <- as.numeric(sample_means[[gi]])
      msel <- grows & partable$matrix == mean_matrix & partable$row %in% x_idx
      if (any(msel)){
        idx <- which(msel)
        partable$est[idx] <- m[partable$row[idx]]
        partable$fixed[idx] <- TRUE
        partable$par[idx] <- 0
      }
    }
  }

  # Renumber the free parameters to be contiguous (1..K):
  freepar <- partable$par != 0
  if (any(freepar)){
    oldpars <- sort(unique(partable$par[freepar]))
    partable$par[freepar] <- match(partable$par[freepar], oldpars)
  }

  partable
}

# Saturated marginal Gaussian log-likelihood of the fixed x-block, summed over
# groups. This is the (constant) likelihood contribution of the exogenous block,
# which is subtracted from the joint log-likelihood so that the reported logl
# matches lavaan's CONDITIONAL log-likelihood under fixed.x (lavaan reports the
# likelihood of the endogenous variables given x). Per group g (sample size n_g,
# ML x-block covariance Sx, p_x variables):
#   ll_x = -n_g/2 ( p_x log(2 pi) + log|Sx| + p_x )
# Uses the stored sample covariances (ML, n denominator) restricted to the x
# variables. Returns 0 when there are no fixed-x variables.
fixed_x_marginal_loglik <- function(x){
  fx <- get_fixed_x(x)
  if (length(fx$idx) == 0) return(0)
  x_idx <- fx$idx
  p_x <- length(x_idx)
  nobs_g <- x@sample@groups$nobs
  total <- 0
  for (g in seq_along(nobs_g)){
    S <- as.matrix(x@sample@covs[[g]])[x_idx, x_idx, drop = FALSE]
    ld <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
    total <- total - nobs_g[g] / 2 * (p_x * log(2 * pi) + ld + p_x)
  }
  total
}

# Fix the exogenous latent block of an lvm parameter table for fixed.x. The named
# exogenous variables are single-indicator latents (loading 1, residual variance
# 0): each maps to exactly one observed indicator. Fixes (per group)
#   - the sigma_zeta block among the exogenous latents to the SAMPLE covariance
#     of their observed indicators (their saturated ML value, since a
#     single-indicator latent's variance equals its indicator's variance), and
#   - the observed indicators' intercepts (nu) to the sample means,
# marking those rows fixed (par = 0) and renumbering the free parameters.
# 'lat_idx' / 'obs_idx' are the 1-based latent / observed indices of the
# exogenous variables (in the same order). Returns the modified parameter table.
apply_fixed_x_partable_lvm <- function(partable, lat_idx, obs_idx,
                                       sample_covs, sample_means, group_ids,
                                       meanstructure){
  if (length(lat_idx) == 0) return(partable)

  for (gi in seq_along(group_ids)){
    gid <- group_ids[gi]
    S <- as.matrix(sample_covs[[gi]])
    grows <- partable$group_id == gid

    # sigma_zeta block among the exogenous latents -> sample cov of indicators:
    for (r in seq_along(lat_idx)){
      for (cc in seq_along(lat_idx)){
        sel <- grows & partable$matrix == "sigma_zeta" &
          partable$row == lat_idx[r] & partable$col == lat_idx[cc]
        if (any(sel)){
          partable$est[sel] <- S[obs_idx[r], obs_idx[cc]]
          partable$fixed[sel] <- TRUE
          partable$par[sel] <- 0
        }
      }
    }

    # Observed exogenous intercepts (nu) -> sample means:
    if (isTRUE(meanstructure) && !is.null(sample_means)){
      m <- as.numeric(sample_means[[gi]])
      for (r in seq_along(obs_idx)){
        sel <- grows & partable$matrix == "nu" & partable$row == obs_idx[r]
        if (any(sel)){
          partable$est[sel] <- m[obs_idx[r]]
          partable$fixed[sel] <- TRUE
          partable$par[sel] <- 0
        }
      }
    }
  }

  freepar <- partable$par != 0
  if (any(freepar)){
    oldpars <- sort(unique(partable$par[freepar]))
    partable$par[freepar] <- match(partable$par[freepar], oldpars)
  }
  partable
}

# Accessor for the stored fixed.x configuration. Returns a list with:
#   names : the exogenous variable names (character(0) if none)
#   idx   : the 1-based indices, into the OBSERVED variable ordering, of the
#           exogenous variables whose marginal block is conditioned on. For
#           varcov these are the named observed variables; for lvm (where
#           fixed_x names single-indicator latents) these are the indices of the
#           corresponding observed indicators, stored at construction time in
#           x@types$fixed_x_obs_idx.
# Read defensively (models created before this feature have no $fixed_x entry).
get_fixed_x <- function(x){
  fx <- x@types$fixed_x
  if (is.null(fx) || length(fx) == 0) return(list(names = character(0), idx = integer(0)))
  obs_idx <- x@types$fixed_x_obs_idx
  if (!is.null(obs_idx) && length(obs_idx) > 0){
    return(list(names = fx, idx = as.integer(obs_idx)))
  }
  # Fallback (varcov): the names are observed variables.
  varlabels <- x@sample@variables$label
  idx <- match(fx, varlabels)
  idx <- idx[!is.na(idx)]
  list(names = fx, idx = idx)
}
