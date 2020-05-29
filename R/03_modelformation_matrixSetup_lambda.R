matrixsetup_lambda <- function(
  lambda, # lambda argument
  nGroup, # Number of groups
  expcov, # Expected covariance matrices (list).
  observednames, latentnames,
  equal = FALSE,
  sampletable,
  name = "lambda",
  identification = c("loadings","variance"),
  simple = FALSE
){
  identification <- match.arg(identification)
  
  # Fix lambda:
  lambda <- fixMatrix(lambda,nGroup,equal)
  
  nLat <- ncol(lambda)
  nObs <- nrow(lambda)
  
  # For each group, form starting values:
  lambdaStart <- lambda
  sigma_epsilon_start <- array(0, dim = c(nObs, nObs, nGroup))
  sigma_zeta_start <- array(0, dim = c(nLat, nLat, nGroup))

  for (g in 1:nGroup){
    # Current cov estimate:
    curcov <- as.matrix(spectralshift(expcov[[g]]))
    if (!(nObs > 3 && nObs > nLat) || simple){
      simple <- TRUE
    } else {
      simple <- FALSE

      tryres <- try({
      
      # Residual and latent varcov:
      fa <- psych::fa(r = curcov, nfactors = nLat, rotate = "promax", covar = TRUE)
      
      # Sigma_epsilon start:
      sigma_epsilon_start[,,g] <- diag(fa$uniquenesses)
      
      load <- fa$loadings
      # Correlation between factor loadings:
      # Generate 100 different permutations and choose the best:
      if (nLat > 1){
        # browser()
        perms <- combinat::permn(1:nLat)
        # fit <- sapply(perms, function(x)sum(diag(stats::cor(abs(load[,x]), lambda[,,g]))))
        fit <- sapply(perms, function(x)sum(abs(load)[,x] * (lambda[,,g]!=0)))
        best <- perms[[which.max(fit)]]        
      } else {
        best <- 1
      }
      
      sigma_zeta_start[,,g] <- fa$r.scores[best,best]
      lambdaStart[,,g] <- lambda[,,g] * load[,best]
      
      # Fix potentially low  and high factor loadings:
      ind <- abs(lambdaStart[,,g]) < 0.25 & abs(lambdaStart[,,g]) > 0
      lambdaStart[,,g][ind] <- sign(lambdaStart[,,g][ind]) * 0.25
      
      ind <- abs(lambdaStart[,,g]) > 2 & abs(lambdaStart[,,g]) > 0
      lambdaStart[,,g][ind] <- sign(lambdaStart[,,g][ind]) * 2
      
      # If loadings identification, I need to rescale ...
      if (identification == "loadings"){
        scaleMat <- matrix(0, nLat, nLat)
        for (i in 1:nLat){
          scaleMat[i,i] <- 1/lambdaStart[,,g,drop=FALSE][,i,1][lambdaStart[,,g,drop=FALSE][,i,1]!=0][1]
        }
        
        lambdaStart[,,g] <- lambdaStart[,,g] %*% scaleMat
        sigma_zeta_start[,,g] <- solve(scaleMat) %*% sigma_zeta_start[,,g] %*% solve(scaleMat)
      }
      })
      if (is (tryres, "try-error")) simple <- TRUE
    }
    
    if (simple){
      sigma_epsilon_start[,,g] <- diag(nObs)
      sigma_zeta_start[,,g] <- diag(nLat)
      lambdaStart[,,g] <- lambda[,,g]
    }

    
# 
#     # For all laodings:
#     for (f in seq_len(ncol(lambdaStart))){
#       # # First principal component of sub cov:
#       # if (any(lambda[,f,g]!=0)){
#       #   ev1 <- eigen(curcov[lambda[,f,g]!=0,lambda[,f,g]!=0])$vectors[,1]
#       #   lambdaStart[lambdaStart[,f,g]!=0,f,g] <- ev1 / ev1[1]        
#       # } 
#       # Univariate factor model:
#       if (any(lambda[,f,g]!=0)){
#         fa <- psych::fa(r = curcov[lambda[,f,g]!=0,lambda[,f,g]!=0], factors = 1, covar = TRUE)
#         load1 <- fa$loadings[,1]
#         
#         if (identification == "loadings"){
#           lambdaStart[lambdaStart[,f,g]!=0,f,g] <- load1 / abs(load1[1])
#         } else {
#           lambdaStart[lambdaStart[,f,g]!=0,f,g] <- load1
#         }
#         
#         # thetaStart[lambdaStart[,f,g]!=0,lambdaStart[,f,g]!=0, g] <- pmin(thetaStart[lambdaStart[,f,g]!=0,lambdaStart[,f,g]!=0, g] , fa$uniquenesses)
#         
#       }
#     }
#     
#  
#     # Now finally:
#     # This means that the factor-part is expected to be:
#     factorPart <- curcov - sigma_epsilon_start[,,g]
#     
#     # Let's take a pseudoinverse:
#     inv <- corpcor::pseudoinverse(kronecker(lambdaStart[,,g],lambdaStart[,,g]))
#     
#     # And obtain psi estimate:
#     sigma_zeta_start[,,g] <- matrix(inv %*% as.vector(factorPart),nLat,nLat)
  }
  
  # Form the model matrix part:
  list(lambda,
       mat = name,
       op =  "~=",
       rownames = observednames,
       colnames = latentnames,
       sparse = TRUE,
       start = lambdaStart,
       sigma_epsilon_start = sigma_epsilon_start,
       sigma_zeta_start = sigma_zeta_start,
       sampletable=sampletable
  )
}