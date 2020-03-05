# Utility functions from Lavaan

# invert positive definite symmetric matrix (eg cov matrix)
# using choleski decomposition
# return log determinant as an attribute
inv.chol <- function(S, logdet=FALSE) {
  cS <- chol(S)
  #if( inherits(cS, "try-error") ) {
  #    print(S)
  #    warning("lavaan WARNING: symmetric matrix is not positive symmetric!")
  #}
  S.inv <- chol2inv( cS )
  
  # Make sparse:
  S.inv[abs(S.inv) < sqrt(.Machine$double.eps)] <- 0
  S.inv <- as(S.inv, "matrix")
  
  if(logdet) {
    diag.cS <- diag(cS)

    # FIXME: Why am I doing - here???
    attr(S.inv, "logdet") <- -sum(log(diag.cS*diag.cS))
  }
  S.inv
}