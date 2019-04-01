expectedmodel <- function(x){
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
  x
}