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
  } else if (x@distribution == "Ising"){
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
