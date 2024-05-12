# This function reduces covariances to limit the correlation to be too strong, useful for starting values (I hope):
maxcor <- function(cov, max = 0.1){
  if (is.list(cov)) {
    return(lapply(cov,maxcor))
  }
  cor <- cov2cor(cov)
  sd <- diag(sqrt(diag(cov)))
  cor[] <- pmin(max,pmax(-max,cor))
  diag(cor) <- 1
  sd %*% cor %*% sd
}
