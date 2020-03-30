# factorstart <- function(covmat, lambda){
#   browser()
#   # Number of latents:
#   nLat <- ncol(lambda)
#   nObs <- ncol(lambda)
#   
#   # Run EFA:
#   fa <- factanal(covmat = covmat, factors = nLat, rotation = "promax")
#   
#   # Uniqueness:
#   thetaStart <- diag(fa$uniquenesses^2)
#   
#   # loadings:
#   loadings <- as.matrix(fa$loadings)
#   class(loadings) <- "matrix"
#   
#   loadings
#   lambda  
#   
#   
# }