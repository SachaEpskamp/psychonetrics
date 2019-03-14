gvar <- function(...){
  var1(...,contemporaneous = "ggm")
}

# # Latent network model creator
# gvar <- function(
#   data, # Dataset
#   beta = "full",
#   omega_zeta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
#   delta_zeta = "full", # If missing, just full for both groups or equal
#   mu,
#   beepvar,
#   dayvar,
#   idvar,
#   vars, # character indicating the variables Extracted if missing from data - group variable
#   groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
#   covs, # alternative covs (array nvar * nvar * ngroup)
#   means, # alternative means (matrix nvar * ngroup)
#   nobs, # Alternative if data is missing (length ngroup)
#   missing = "pairwise",
#   equal = "none", # Can also be any of the matrices
#   baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
#   # fitfunctions, # Leave empty
#   estimator = "ML",
#   optimizer = "default", #ucminf",
#   rawts = FALSE
# ){
#   
#   # FIXME: Not sure why needed...
#   if (missing(vars)) vars2 <- NULL else vars2 <- vars
#   if (missing(idvar)) idvar <- NULL
#   if (missing(dayvar)) dayvar <- NULL
#   if (missing(beepvar)) beepvar <- NULL
#   if (missing(groups)) groups <- NULL
#   
#   # If data is missing with rawts, stop:
#   if (rawts && missing(data)){
#     stop("'data' may not be missing with rawts = TRUE")
#   }
#   if (!missing(data)){
#     data <- as.data.frame(data)
#     if (is.null(names(data))){
#       names(data) <- paste0("V",seq_len(ncol(data)))
#     }
#   }
#   
#   # If data is not missing, make augmented data:
#   if (!missing(data) && !rawts){
#     data <- tsData(data, vars = vars2, beepvar = beepvar, dayvar = dayvar, idvar = idvar, groupvar = groups)
#   }
#   # Extract var names:
#   if (is.null(groups)){
#     vars <- colnames(data)  
#   } else {
#     vars <- colnames(data)[colnames(data)!=groups]
#   }
#   
#   # Obtain sample stats:
#   sampleStats <- samplestats(data = data, 
#                              vars = vars, 
#                              groups = groups,
#                              covs = covs, 
#                              means = means, 
#                              nobs = nobs, 
#                              missing  = ifelse(estimator == "FIML","pairwise",missing),
#                              rawts = rawts,
#                              fimldata = estimator == "FIML")
#   
#   # Check if number of variables is an even number:
#   if (!rawts && nrow(sampleStats@variables)%%2!=0){
#     stop("Number of variables is not an even number: variance-covariance matrix cannot be a Toeplitz matrix. ")
#   }
#   
#   # Check some things:
#   if (rawts){
#     nNode <- nrow(sampleStats@variables) 
#   } else {
#     nNode <- nrow(sampleStats@variables) / 2
#   }
#   
#   # Generate model object:
#   model <- generate_psychonetrics(model = "gvar",sample = sampleStats,computed = FALSE, 
#                                   equal = equal,rawts = rawts,
#                                   optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
#   
#   # Number of groups:
#   nGroup <- nrow(model@sample@groups)
#   
#   # FIXME: Keep this the same for now for rawts = TRUE
#   nVar <- nNode * 2
#   # Add number of observations:
#   model@sample@nobs <-  
#     nVar * (nVar+1) / 2 * nGroup + # Covariances per group
#     nVar * nGroup # Means per group
#   
#   # Fix mu
#   if (rawts){
#     mu <- fixMu(mu,nGroup,nNode,"mu" %in% equal)  
#   } else {
#     mu <- fixMu(mu,nGroup,nNode*2,"mu" %in% equal)
#   }
#   
#   
#   # Fix the omega_zeta matrix:
#   omega_zeta <- fixAdj(omega_zeta,nGroup,nNode,"omega_zeta" %in% equal,diag0=TRUE)
#   
#   # Fix delta:
#   delta_zeta <- fixAdj(delta_zeta,nGroup,nNode,"delta_zeta" %in% equal,diagonal=TRUE)
#   
#   # Fix beta matrix:
#   beta <- fixMatrix(beta,nGroup,nNode,nNode,"beta" %in% equal)
#   
#   # Exogenous varcov:
#   if (!rawts){
#     exoVarCov <- fixAdj("full",nGroup,nNode)    
#   }
#   
#   
#   # Generate starting values:
#   omegaStart <- omega_zeta
#   deltaStart <- delta_zeta
#   muStart <- mu
#   betaStart <- beta
#   
#   if (!rawts){
#     exoStart <- exoVarCov    
#   }
#   
#   
#   # If we are in the rawts world.. let's get la cov matrix estimate from lavaan:
#   if (rawts){
#     tempdata <- tsData(data, vars = vars2, beepvar = beepvar, dayvar = dayvar, idvar = idvar, groupvar = groups)
#     if (missing(groups) || is.null(groups)){
#       lavOut <- lavCor(as.data.frame(tempdata), missing = "fiml", output="lavaan")
#       lavOut <- lavaan::lavInspect(lavOut, what = "sample")
#       curmeans <- list(lavOut$mean)
#       curcovs <- list(lavOut$cov)
#     } else {
#       lavOut <- lavCor(as.data.frame(tempdata), missing = "fiml", output="lavaan", group = group)
#       lavOut <- lavaan::lavInspect(lavOut, what = "sample")
#       curmeans <- lavOut$mean
#       curcovs <- lavOut$cov
#     }
#     
#   } else {
#     curmeans <- sampleStats@means
#     curcovs <- sampleStats@covs
#   }
#   
#   for (g in 1:nGroup){
#     
#     # If rawts, make some dummy means and varcovs:
#     curcovs[[g]] <- as.matrix(spectralshift(curcovs[[g]]))
#     curSstar <-  curcovs[[g]][1:nNode,1:nNode]
#     curMeans <-  curmeans[[g]]
#     curS0 <- curcovs[[g]][nNode + (1:nNode),nNode + (1:nNode)]
#     curS1 <- curcovs[[g]][nNode + (1:nNode),1:nNode]
#     
# 
#     # Means with sample means:
#     if (rawts){
#       muStart[,g] <- curMeans[nNode + (1:nNode)]
#     } else {
#       muStart[,g] <- curMeans    
#     }
#     
#     
#     # Let's get three blocks:
#     if (!rawts){
#       # exoStart[,,g] <- curSstar  
#       # Cholesky decomposition
#       Lest <- t(as.matrix(chol(curSstar)))
#       exoStart[,,g] <- Lest    
#     }
#     
#     S1 <- curS1
#     S0 <- curS0
#     S0inv <- corpcor::pseudoinverse(S0)
#     
#     # A prior guess for beta is:
#     betaStart[,,g] <- (beta[,,g]!=0) * as.matrix(S1 %*% S0inv)
#     
#     # A prior guess for the contemporaneous covariances is (Schur complement):
#     contCov <- as.matrix(curSstar - t(S1) %*% S0inv %*% S1)
#     
#     # Make posdef if not posdev:
#     if (any(eigen(contCov)$values < 0)){
#       contCov <- as.matrix(Matrix::nearPD(contCov)$mat)
#     }
#     
#     # Let's use this as starting estimate:
#     zeroes <- which(omega_zeta[,,g]==0 & t(omega_zeta[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
#     
#     
#     
#     if (nrow(zeroes) == 0){
#       wi <- glasso(contCov, rho = 0.1)$wi
#     } else {
#       glas <- glasso(contCov,
#                      rho = 0.1, zero = zeroes)
#       wi <- glas$wi
#     }
#     wi <- 0.5*(wi + t(wi))
#     wi <- spectralshift(wi)
#     
#     # Network starting values:
#     pcors <- as.matrix(qgraph::wi2net(as.matrix(wi)))
#     pcors[upper.tri(pcors)] <- 0
#     omegaStart[,,g] <- ifelse(pcors!=0,pcors,ifelse(omegaStart[,,g]!=0,0.1,0))
#     # omegaStart[,,g] <- 1*(omegaStart[,,g]!=0) * 0.05
#     diag(omegaStart[,,g] ) <- 0
#     # Delta:
#     # deltaStart[,,g] <- diag(1/sqrt(diag(wi)))
#     deltaStart[,,g] <- diag(sqrt(diag(wi)))
#     
#     # If we are in the FIML world, all this stuff is dangerous, and I use conservative starting values instead:
#     # if (estimator == "FIML"){
#     #   deltaStart[,,g] <- ifelse(deltaStart[,,g]!=0,1,0)
#     #   omegaStart[,,g] <- 1*(omegaStart[,,g]!=0) * ifelse(omegaStart[,,g]>0,0.01,-0.01)
#     #   betaStart[,,g] <- 1*(betaStart[,,g]!=0) * ifelse(betaStart[,,g]>0,0.01,-0.01)
#     #   
#     #   if (!rawts){
#     #     exoStart[,,g] <- 1*(exoStart[,,g]!=0) * ifelse(exoStart[,,g]>0,0.01,-0.01)
#     #     diag(exoStart[,,g]) <- 1
#     #   }
#     #   
#     # }
#   }
#   
#   # Full parameter table:
#   matList <- list()
#   matList$mu <-  list(mu,
#                       mat =  "mu",
#                       op =  "~1",
#                       symmetrical= FALSE, 
#                       sampletable=sampleStats,
#                       rownames = sampleStats@variables$label,
#                       colnames = "1",
#                       start = muStart)
#   
#   if (!rawts){
#     matList$exo_cholesky <-  list(exoVarCov,
#                                mat =  "exo_cholesky",
#                                op =  "~chol~",
#                                lowertri = TRUE, 
#                                sampletable=sampleStats,
#                                rownames = sampleStats@variables$label[1:nNode],
#                                colnames = sampleStats@variables$label[1:nNode],
#                                start = exoStart)    
#   }
#   
#   matList$beta <-  list(beta,
#                         mat =  "beta",
#                         op =  "<-",
#                         sampletable=sampleStats,
#                         rownames = sampleStats@variables$label[((!rawts)*nNode) + (1:nNode)],
#                         colnames = sampleStats@variables$label[1:nNode],
#                         sparse = TRUE,
#                         start = betaStart
#   )
#   
#   matList$omega_zeta <- list(omega_zeta,
#                              mat =  "omega_zeta",
#                              op =  "--",
#                              symmetrical= TRUE, 
#                              sampletable=sampleStats,
#                              rownames = sampleStats@variables$label[((!rawts)*nNode) + (1:nNode)],
#                              colnames = sampleStats@variables$label[((!rawts)*nNode) + (1:nNode)],
#                              sparse = TRUE,
#                              posdef = TRUE,
#                              diag0=TRUE,
#                              lower = -1,
#                              upper = 1,
#                              start = omegaStart
#   )
#   
#   matList$delta_zeta <-  list(delta_zeta,
#                               mat =  "delta_zeta",
#                               op =  "~/~",
#                               symmetrical= TRUE, 
#                               sampletable=sampleStats,
#                               rownames = sampleStats@variables$label[((!rawts)*nNode) + (1:nNode)],
#                               colnames = sampleStats@variables$label[((!rawts)*nNode) + (1:nNode)],
#                               sparse = TRUE,
#                               posdef = TRUE,
#                               diagonal = TRUE,
#                               lower = 0.01,
#                               start = deltaStart
#   )
#   
#   # Generate the full parameter table:
#   pars <- do.call(generateAllParameterTables, matList)
#   
#   
#   # Store in model:
#   model@parameters <- pars$partable
#   model@matrices <- pars$mattable
#   model@extramatrices <- list(
#     D =  psychonetrics::duplicationMatrix(nNode*2), # Toeplitz matrix D 
#     D2 = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
#     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
#     Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
#     In = Diagonal(nNode), # Identity of dim n
#     In2 = Diagonal(nNode), # Identity of dim n^2
#     A = psychonetrics::diagonalizationMatrix(nNode),
#     C = as(lavaan::lav_matrix_commutation(nNode,nNode),"sparseMatrix")
#     # P=P # Permutation matrix
#   )
#   
#   
#   if (!rawts){
#     # Come up with P...
#     # Dummy matrix to contain indices:
#     dummySigma <- matrix(0,nNode*2,nNode*2)
#     smallMat <- matrix(0,nNode,nNode)
#     dummySigma[1:nNode,1:nNode][lower.tri(smallMat,diag=TRUE)] <- seq_len(nNode*(nNode+1)/2)
#     dummySigma[nNode + (1:nNode),nNode + (1:nNode)][lower.tri(smallMat,diag=TRUE)] <- max(dummySigma) + seq_len(nNode*(nNode+1)/2)
#     dummySigma[nNode + (1:nNode),1:nNode] <- max(dummySigma) + seq_len(nNode^2)
#     inds <- dummySigma[lower.tri(dummySigma,diag=TRUE)]
#     
#     # P matrix:
#     # P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=inds))
#     model@extramatrices$P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=order(inds)))
#   }
#   
#   
#   
#   
#   
#   # If we are in the rawts world, we need Drawts duplication matrix for every group...
#   if (rawts){
#     # I need to fill this list:
#     Drawts <- list()
#     
#     # Dummy matrices with indices:
#     muDummy <- matrix(1:nNode)
#     sigDummy <- matrix(0,nNode,nNode)
#     sigDummy[lower.tri(sigDummy,diag=TRUE)] <- max(muDummy) + seq_len(nNode*(nNode+1)/2)
#     sigDummy[upper.tri(sigDummy)] <- t(sigDummy)[upper.tri(sigDummy)]
#     
#     U <- list(sigDummy)
#     # Now make all lag-k blocks...
#     maxN <- max(sapply(sampleStats@missingness,nrow))
#     
#     # Form blocks:
#     for (i in 1:(maxN-1)){
#       U[[length(U) + 1]] <- matrix(max(unlist(U)) + seq_len(nNode^2), nNode, nNode)
#     }
#     
#     # # Magic code from stackovervlow:
#     #   k <- min(unlist(lapply(U, dim)))
#     #   n <- length(U)
#     #   #
#     #   # Create the "strip".
#     #   #
#     #   strip <- array(NA, dim=c(k,k,2*n-1))
#     #   for (i in 1:n) strip[,,i] <- U[[n+1-i]]
#     #   if (n > 1) for (i in 2:n) strip[,,n+i-1] <- t(U[[i]])
#     #   #
#     #   # Assemble into "block-Toeplitz" form.
#     #   #
#     #   allSigmas <- array(NA, dim=c(k,k,n,n))
#     #   #
#     #   # Blast the strip across X.
#     #   #
#     #   for (i in 1:n) allSigmas[,,,i] <- strip[,,(n+1-i):(2*n-i)]
#     #   allSigmas <- matrix(aperm(allSigmas, c(1,3,2,4)), n*k)
#     allSigmas <- blockToeplitz(U)
#     
#     
#     for (g in 1:nGroup){
#       missings <- sampleStats@missingness[[g]]
#       
#       # Create the massive matrix:
#       muFull <- Reduce("rbind",lapply(seq_len(nrow(missings)),function(x)muDummy))
#       
#       # sigFull <- Reduce("bdiag",lapply(seq_len(nrow(missings)),function(x)sigDummy))
#       sigFull <- allSigmas[seq_len(nrow(missings)*nNode), seq_len(nrow(missings)*nNode)]
#       
#       # Cut out the rows and cols:
#       obsvec <- !as.vector(t(missings))
#       muFull <- muFull[obsvec,]
#       sigFull <- as.matrix(sigFull[obsvec,obsvec])
#       
#       # Now make a vector with the indices for the full rawts distiributional parameter set:
#       distVec <- c(as.vector(muFull),as.vector(sigFull))
#       nTotal <- length(distVec)
#       distVecrawts <- seq_along(distVec)[distVec!=0]
#       distVec <- distVec[distVec!=0]
#       # Now I can make the matrix:
#       
#       Drawts[[g]] <- sparseMatrix(
#         i = distVecrawts, j = distVec, dims = c(nTotal, max(U[[nrow(missings)]]))
#       )
#       
#     }
#     
#     
#     # Add these to the model:
#     model@Drawts = Drawts
#   }
#   
#   # Form the model matrices
#   model@modelmatrices <- formModelMatrices(model)
#   
#   
#   ### Baseline model ###
#   if (baseline_saturated){
#     
#     
#     # Normally via ggm:
#     if (!rawts){
#       # Baseline GGM should be block matrix:
#       basGGM <- diag(nNode*2)
#       basGGM[1:nNode,1:nNode] <- 1
#       
#       model@baseline_saturated$baseline <- cholesky(data = data,
#                                                     lowertri = basGGM,
#                                                     vars = vars,
#                                                     groups = groups,
#                                                     covs = covs,
#                                                     means = means,
#                                                     nobs = nobs,
#                                                     missing = missing,
#                                                     equal = equal,
#                                                     estimator = estimator,
#                                                     baseline_saturated = FALSE)
#       
#       # Add model:
#       # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
#       ### Saturated model ###
#       model@baseline_saturated$saturated <- cholesky(data = data, 
#                                                      lowertri = "full", 
#                                                      vars = vars,
#                                                      groups = groups,
#                                                      covs = covs,
#                                                      means = means,
#                                                      nobs = nobs,
#                                                      missing = missing,
#                                                      equal = equal,
#                                                      estimator = estimator,
#                                                      baseline_saturated = FALSE)
#       
#       # if not FIML, Treat as computed:
#       if (estimator != "FIML"){
#         model@baseline_saturated$saturated@computed <- TRUE
#         
#         # FIXME: TODO
#         model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)      
#       }
#     } else {
#       
#       # Via gvar:
#       model@baseline_saturated$baseline <- gvar(data,omega_zeta = "empty", beta = "empty",
#                                                 vars = vars,
#                                                 groups = groups,
#                                                 missing = missing,
#                                                 equal = equal,
#                                                 baseline_saturated = FALSE, rawts = TRUE)
#       
#       # Add model:
#       # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
#       
#       
#       ### Saturated model ###
#       model@baseline_saturated$saturated <- gvar(data,omega_zeta = "full", beta = "full",
#                                                  vars = vars,
#                                                  groups = groups,
#                                                  missing = missing,
#                                                  equal = equal,
#                                                  baseline_saturated = FALSE, rawts = TRUE)
#       
#       
#       # Treat as computed:
#       model@baseline_saturated$saturated@computed <- TRUE
#       
#       # Add saturated fit
#       model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
#       
#     }
#     
#   }
#   
#   
#   # Return model:
#   return(model)
# }
