start <- function(mat,list,alt){
  if (!is.null(list[[mat]])){
    return(list[[mat]])
  } else {
    return(alt)
  }
}

# toFree: function to take input and return NA whenever NA or character:
toFree <- function(x){
  suppressWarnings(mode(x) <- "numeric")
  x
}
isFree <- function(x) is.na(toFree(x))
# to label, convert to equality constraints:
toLabel <- function(x, mat, symmetric=FALSE, group = 1){
  isNA <- is.na(x)
  free <- toFree(x)
  suppressWarnings(mode(x) <- "character")
  inds <- which(isNA | !is.na(free),arr.ind=TRUE)
  x[isNA | !is.na(free)] <- paste0(mat,"_",inds[,1],"_",inds[,2],"_",group)
  if (symmetric){
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
  }
  # x[isNA | !is.na(free)]  <- NA
  return(x)
}

# And from graphicalVAR:
computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(Matrix::forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  if (ncol(beta) == nrow(beta)+1){
    beta <- beta[,-1,drop=FALSE]
  }
  sigma <- solve(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}
