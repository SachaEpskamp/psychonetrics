blockToeplitz <- function(U){
  # Magic code from stackovervlow:
  k <- min(unlist(lapply(U, dim)))
  n <- length(U)
  #
  # Create the "strip".
  #
  strip <- array(NA, dim=c(k,k,2*n-1))
  for (i in 1:n) strip[,,i] <- t(U[[n+1-i]])
  if (n > 1) for (i in 2:n) strip[,,n+i-1] <- U[[i]]
  #
  # Assemble into "block-Toeplitz" form.
  #
  X <- array(NA, dim=c(k,k,n,n))
  #
  # Blast the strip across X.
  #
  for (i in 1:n) X[,,,i] <- strip[,,(n+1-i):(2*n-i)]
  X <- matrix(aperm(X, c(1,3,2,4)), n*k)
  X
}