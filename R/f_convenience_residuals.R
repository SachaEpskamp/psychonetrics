psychonetrics_residuals <- function(object, round = TRUE, digits = 2){
  obsCovs <- object@sample@covs
  obsCors <- object@sample@cors
  obsMeans <- object@sample@means
  
  impCovs <- lapply(object@modelmatrices,"[[","sigma")
  impCors <- lapply(impCovs, function(x)cov2cor(as.matrix(x)))
  impMeans <-   lapply(object@modelmatrices,"[[","mu")
  
  # Final list:
  nGroup <- nrow(object@sample@groups)
  res <- vector("list", nGroup)
  for (i in seq_len(nGroup)){
    res[[i]] <- list(
      covs = as.matrix(obsCovs[[i]] - impCovs[[i]]),
      cors = as.matrix(obsCors[[i]] - impCors[[i]]),
      means = as.matrix(obsMeans[[i]] - impMeans[[i]])
    )
    if (round){
      res[[i]] <- lapply(res[[i]],round,digits=digits)
    }
  }
  names(res) <- object@sample@groups$label
  return(res)
}

setMethod(f = "resid",
          signature = "psychonetrics",
          definition = psychonetrics_residuals)

setMethod(f = "residuals",
          signature = "psychonetrics",
          definition = psychonetrics_residuals)