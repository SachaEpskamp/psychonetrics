expectedmodel <- function(x){
  if (x@distribution == "Gaussian"){
    start <- parVector(x)
    prep <- prepareModel(start, x)
    for (g in 1:nrow(x@sample@groups)){
      x@sample@means[[g]] <- prep$groupModels[[g]]$mu
      x@sample@covs[[g]] <- prep$groupModels[[g]]$sigma

      if (length(x@sample@fimldata) > 0){
        nPat <- length(x@sample@fimldata[[g]])
        for (i in seq_len(nPat)){
          x@sample@fimldata[[g]][[i]]$means <- as.matrix(prep$groupModels[[g]]$mu[!x@sample@fimldata[[g]][[i]]$pattern,drop=FALSE])
          if (!all(x@sample@fimldata[[g]][[i]]$S == 0)){
            x@sample@fimldata[[g]][[i]]$S <- as.matrix(prep$groupModels[[g]]$sigma[!x@sample@fimldata[[g]][[i]]$pattern,!x@sample@fimldata[[g]][[i]]$pattern, drop = FALSE])
          }
        }
      }
    }
  } else if (x@distribution == "TwoLevelGaussian"){
    # Two-level ML with WITHIN-cluster MISSING data: the sufficient statistics
    # do not compress and there is no per-size moment to replace by its
    # expectation. numeric_FisherInformation() then differentiates the
    # gradient at the given parameters on the raw data, i.e. it uses the
    # (numeric) observed-information-based unit information -- a well-defined
    # estimator of the information. So leave the model unchanged here:
    if (twolevel_model_has_missing(x)){
      return(x)
    }
    # Two-level ML estimator (ml_lvm), COMPLETE data: replace the two-level
    # sufficient statistics by their model-implied expectations, so that the
    # analytic gradient is zero in expectation and its numeric Jacobian equals
    # the expected (Fisher) information. With Omega_s = Sigma_B + Sigma_W / n_s
    # the expectations are E(S_PW) = Sigma_W, E(mean.d_s) = mu and
    # E(A_s) = Omega_s (split here as cov.d_s = Omega_s, mean.d_s = mu):
    start <- parVector(x)
    prep <- prepareModel(start, x)
    twolevelStats <- get_twolevel_stats(x@sample)
    for (g in 1:nrow(x@sample@groups)){
      mu_g <- as.vector(prep$groupModels[[g]]$mu)
      SW_g <- as.matrix(prep$groupModels[[g]]$sigma_within)
      SB_g <- as.matrix(prep$groupModels[[g]]$sigma_between)
      twolevelStats[[g]]$S_PW <- SW_g
      sizes <- twolevelStats[[g]]$sizes
      for (s in seq_len(nrow(sizes))){
        twolevelStats[[g]]$mean_d[[s]] <- mu_g
        twolevelStats[[g]]$cov_d[[s]] <- SB_g + SW_g / sizes$nj[s]
      }
    }
    x@sample@twolevel <- twolevelStats
  } else if (x@distribution == "Spin"){
    start <- parVector(x)
    prep <- prepareModel(start, x)
    for (g in 1:nrow(x@sample@groups)){
      nobs_g <- x@sample@groups$nobs[g]
      # Replace sample means with model-implied E(y_i):
      x@sample@means[[g]] <- as.matrix(prep$groupModels[[g]]$exp_v1)
      # Replace sample squares (sum of products) with model-implied nobs * E(y_i * y_j):
      x@sample@squares[[g]] <- prep$groupModels[[g]]$exp_v2 * nobs_g
    }
  }


  return(x)
}
