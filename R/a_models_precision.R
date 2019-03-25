# Latent network model creator
precision <- prec <- function(
  ...
){
 model <- varcov(...,type = "prec")
  
  # Return model:
  return(model)
}
# # Latent network model creator
# ggm <- function(
#   data, # Dataset
#   omega = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
#   delta = "full", # If missing, just full for both groups or equal
#   mu,
#   vars, # character indicating the variables Extracted if missing from data - group variable
#   groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
#   covs, # alternative covs (array nvar * nvar * ngroup)
#   means, # alternative means (matrix nvar * ngroup)
#   nobs, # Alternative if data is missing (length ngroup)
#   missing = "listwise",
#   equal = "none", # Can also be any of the matrices
#   baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
#   estimator = "ML",
#   optimizer = "default", #"ucminf",
#   rawts = FALSE # Set to TRUE to do rawts estimation instead...
# ){
#   if (rawts){
#     warning("'rawts' is only included for testing purposes. Please do not use!")
#   }
# 
#   # Obtain sample stats:
#   sampleStats <- samplestats(data = data, 
#                              vars = vars, 
#                              groups = groups,
#                              covs = covs, 
#                              means = means, 
#                              nobs = nobs, 
#                              missing = ifelse(estimator == "FIML","pairwise",missing),
#                              rawts = rawts,
#                              fimldata = estimator == "FIML")
# 
#   # Check some things:
#   nNode <- nrow(sampleStats@variables)
#   
#   # Generate model object:
#   model <- generate_psychonetrics(model = "ggm",sample = sampleStats,computed = FALSE, 
#                                   equal = equal,
#                                   optimizer = optimizer, estimator = estimator, distribution = "Gaussian",
#                                   rawts = rawts)
#   
#   # Number of groups:
#   nGroup <- nrow(model@sample@groups)
#   
#   # Add number of observations:
#   model@sample@nobs <-  
#     nNode * (nNode+1) / 2 * nGroup + # Covariances per group
#     nNode * nGroup # Means per group
#   
#   # Fix mu
#   mu <- fixMu(mu,nGroup,nNode,"mu" %in% equal)
#   
#   # Fix omega
#   omega <- fixAdj(omega,nGroup,nNode,"omega" %in% equal,diag0=TRUE)
#   
#   # Fix delta:
#   delta <- fixAdj(delta,nGroup,nNode,"delta" %in% equal,diagonal=TRUE)
#   
#   
#   # Generate starting values:
#   omegaStart <- omega
#   deltaStart <- delta
#   muStart <- mu
#   
#   # If we are in the rawts world.. let's get la cov matrix estimate from lavaan:
#   if (rawts){
#     if (missing(groups)){
#       lavOut <- lavCor(as.data.frame(data[,vars]), missing = "fiml", output="lavaan")
#     } else {
#       lavOut <- lavCor(as.data.frame(data[,vars]), missing = "fiml", output="lavaan", group = group)
#     }
#     lavOut <- lavaan::lavInspect(lavOut, what = "sample")
#   }
#   
#   for (g in 1:nGroup){
#     
#     # Are we in the rawts world?
#     if (rawts){
#       if (missing(groups)){
#         covest <-as.matrix( lavOut$cov)
#         meanest <- lavOut$mean
#       } else {
#         covest <- as.matrix(lavOut$cov[[g]])
#         meanest <- lavOut$mean[[g]]
#       }
#     } else {
#       covest <- as.matrix(sampleStats@covs[[g]])
#       meanest <-  sampleStats@means[[g]]
#     }
#       
#       # Means with sample means:
#       muStart[,g] <- ifelse(is.na(meanest),0,meanest)
#       
#       # If the estimator is fiml, let's not bother with the sample var cov matrix and instead start with some general starting values:
#       # if (estimator == "fiml"){
#       #   omegaStart[,,g] <- 1*(omegaStart[,,g]!=0) * ifelse(covest > 0, 0.05, -0.05)
#       #   deltaStart[,,g] <- 1*(deltaStart[,,g]!=0) * diag(nrow( deltaStart[,,g]))
#       # } else {
#         # For the network, compute a rough glasso network:
#         zeroes <- which(omegaStart[,,g]==0 & t(omegaStart[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
#         if (nrow(zeroes) == 0){
#           wi <- solve_symmetric(covest)
#         } else {
#           glas <- glasso(as.matrix(covest),
#                          rho = 1e-10, zero = zeroes)
#           wi <- glas$wi
#         }
#         
#         # Network starting values:
#         omegaStart[,,g] <- as.matrix(qgraph::wi2net(wi))
#         diag(omegaStart[,,g] ) <- 0
#         
#         # # FIXME: Quick and dirty fiml fix:
#         # if (estimator == "FIML"){
#         #   omegaStart[,,g] <- 0.5 * omegaStart[,,g]
#         # }
#         
#         # Delta:
#         deltaStart[,,g] <- diag(1/sqrt(diag(wi)))
#       # }
# 
#   }
# 
#   # Generate the full parameter table:
#   pars <- generateAllParameterTables(
#     # Mu:
#     list(mu,
#          mat =  "mu",
#          op =  "~1",
#          symmetrical= FALSE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = "1",
#          start = muStart),
#     
#  
#     # Omega:
#     list(omega,
#          mat =  "omega",
#          op =  "--",
#          symmetrical= TRUE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = sampleStats@variables$label,
#          sparse = TRUE,
#          posdef = TRUE,
#          diag0=TRUE,
#          lower = -1,
#          upper = 1,
#          start = omegaStart
#     ),
#     
#     # Delta_eta:
#     list(delta,
#          mat =  "delta",
#          op =  "~/~",
#          symmetrical= TRUE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = sampleStats@variables$label,
#          sparse = TRUE,
#          posdef = TRUE,
#          diagonal = TRUE,
#          lower = 0,
#          start = deltaStart
#     )
#     
#   )
#   
#   # Store in model:
#   model@parameters <- pars$partable
#   model@matrices <- pars$mattable
#   model@extramatrices <- list(
#     D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
#     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
#     Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
#     In = Diagonal(nNode), # Identity of dim n
#     A = psychonetrics::diagonalizationMatrix(nNode)
#   )
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
#     muDummy <- as(muDummy,"Matrix")
#     sigDummy <- as(sigDummy,"Matrix")
#     for (g in 1:nGroup){
#       missings <- sampleStats@missingness[[g]]
#       
#       # Create the massive matrix:
#       muFull <- Reduce("rbind",lapply(seq_len(nrow(missings)),function(x)muDummy))
#       sigFull <- Reduce("bdiag",lapply(seq_len(nrow(missings)),function(x)sigDummy))
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
#         i = distVecrawts, j = distVec, dims = c(nTotal, max(sigDummy))
#       )
#     }
# 
# 
#     # Add these to the model:
#     model@Drawts = Drawts
#   }
#   
#   
#   # Form the model matrices
#   model@modelmatrices <- formModelMatrices(model)
#   
#   
#   ### Baseline model ###
#   if (baseline_saturated){
#     
#     # Form baseline model:
#     model@baseline_saturated$baseline <- cholesky(data, 
#                                                 lowertri = "empty",
#                                                 vars = vars,
#                                                 groups = groups,
#                                                 covs = covs,
#                                                 means = means,
#                                                 nobs = nobs,
#                                                 missing = missing,
#                                                 equal = equal,
#                                                 estimator = estimator,
#                                                 baseline_saturated = FALSE)
#     
#     # Add model:
#     # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
#     
#     
#     ### Saturated model ###
#     model@baseline_saturated$saturated <- cholesky(data, 
#                                                  lowertri = "full", 
#                                                  vars = vars,
#                                                  groups = groups,
#                                                  covs = covs,
#                                                  means = means,
#                                                  nobs = nobs,
#                                                  missing = missing,
#                                                  equal = equal,
#                                                  estimator = estimator,
#                                                  baseline_saturated = FALSE)
#     
#     # if not FIML, Treat as computed:
#     if (estimator != "FIML"){
#       model@baseline_saturated$saturated@computed <- TRUE
#       
#       # FIXME: TODO
#       model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)      
#     }
#   }
#   
#   # Return model:
#   return(model)
# }