sem <- function(...){
  lvm(...,latent = "cov", residual = "cov")
}

# # Latent network model creator
# rnm <- function(
#   data, # Dataset
#   lambda, # only required non-missing matrix
#   sigma_zeta = "auto", # full for exogenous, diagonal for endogenous
#   omega_epsilon = "empty", # Residuals, will default to a diagonal matrix
#   delta_epsilon = "full",
#   beta = "empty",
#   tau,
#   tau_eta,
#   vars, # character indicating the variables Extracted if missing from data - group variable
#   latents, # Name of latent varianles
#   groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
#   covs, # alternative covs (array nvar * nvar * ngroup)
#   means, # alternative means (matrix nvar * ngroup)
#   nobs, # Alternative if data is missing (length ngroup)
#   missing = "listwise",
#   equal = "none", # Can also be any of the matrices
#   baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
#   # fitfunctions, # Leave empty
#   identify = TRUE,
#   identification = c("loadings","variance"),
#   estimator = "ML",
#   optimizer = "default"
# ){
#   identification <- match.arg(identification)
#   
#   # Obtain sample stats:
#   sampleStats <- samplestats(data = data, 
#                              vars = vars, 
#                              groups = groups,
#                              covs = covs, 
#                              means = means, 
#                              nobs = nobs, 
#                              missing  = ifelse(estimator == "FIML","pairwise",missing),
#                              fimldata = estimator == "FIML")
#   
#   # Check some things:
#   nNode <- nrow(sampleStats@variables)
#   
#   # Generate model object:
#   model <- generate_psychonetrics(model = "rnm",sample = sampleStats,computed = FALSE, 
#                                   equal = equal,identification = identification,
#                                   optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
#   
#   # Number of groups:
#   nGroup <- nrow(model@sample@groups)
#   
#   # Add number of observations:
#   model@sample@nobs <-  
#     nNode * (nNode+1) / 2 * nGroup + # Covariances per group
#     nNode * nGroup # Means per group
#   
#   # Stop if lambda is missing or a character:
#   if (missing(lambda)){
#     stop("'lambda' may not be missing")
#   }
#   if (is.character(lambda)){
#     stop("'lambda' may not be a string")
#   }
#   
#   # Fix tau
#   tau <- fixMu(tau,nGroup,nNode,"tau" %in% equal)
#   
#   # Fix the lambda matrix:
#   lambda <- fixMatrix(lambda,nGroup,equal = "lambda" %in% equal)
#   
#   # Number of latents:
#   nLatent <- dim(lambda)[2]
#   
#   # Fix tau_eta
#   tau_eta <- fixMu(tau_eta,nGroup,nLatent,"tau_eta" %in% equal)
#   
#   # Fix the beta matrix:
#   beta <- fixMatrix(beta,nGroup,nLatent,nLatent,"beta" %in% equal,diag0=TRUE)
# 
#   # Fix sigma_zeta:
#   if (is.character(sigma_zeta) && sigma_zeta == "empty"){
#     diag0 <- TRUE
#   } else {
#     diag0 <- FALSE
#   }
#   if (is.character(sigma_zeta) && sigma_zeta == "auto"){
#     autoSigma_zeta <- TRUE
#     sigma_zeta <- "diag"
#   } else {
#     autoSigma_zeta <- FALSE
#   }
#   if (is.character(sigma_zeta) && sigma_zeta == "diag"){
#     diagonal <- TRUE
#   } else {
#     diagonal <- FALSE
#   }
#   sigma_zeta <- fixAdj(sigma_zeta,nGroup,nLatent,"sigma_zeta" %in% equal,diag0 = diag0,diagonal = diagonal)
#   # Auto stuff:
#   if (autoSigma_zeta){
#     # for every group:
#     for (g in 1:nGroup){
#       # Which are exogenous?
#       exo <- rowSums(beta[,,g]!=0)==0
#       sigma_zeta[exo,exo,g] <- ifelse(lower.tri(sigma_zeta[exo,exo,g],diag=TRUE),1,0)
#     }
#     
#     # Reset eqality constrains if needed:
#     if ("sigma_zeta" %in% equal){
#       for (i in 1:nrow(sigma_zeta)){
#         for (j in 1:i){
#           if (sigma_zeta[i,j,1]!=0){
#             sigma_zeta[i,j,] <- max(sigma_zeta) + 1
#           }
#         }
#       }
#     }
#   }
#   
#   # Fix delta_epsilon
#   if (is.character(delta_epsilon) && delta_epsilon == "empty"){
#     diag0 <- TRUE
#   } else {
#     diag0 <- FALSE
#   }
#   delta_epsilon <- fixAdj(delta_epsilon,nGroup,nNode,"delta_epsilon" %in% equal,diagonal=TRUE,diag0=diag0)
#   
#   # Fix omega_epsilon:
#   omega_epsilon <- fixAdj(omega_epsilon,nGroup,nNode,"omega_epsilon" %in% equal,diag0=TRUE)
#   
#   # If latents is not provided, make it:
#   if (missing(latents)){
#     latents <- paste0("Eta_",seq_len(nLatent))
#   }
#   if (length(latents) != nLatent){
#     stop("Length of 'latents' is not equal to number of latent variables in model.")
#   }
#   
#   
#   # Generate starting values:
#   lambdaStart <- lambda
#   tauStart <- tau
#   sigma_zeta_start <- sigma_zeta
#   delta_epsilon_start <- delta_epsilon
#   betaStart <- 0.1*(beta!=0)
#   for (g in 1:nGroup){
#     # Means with sample means:
#     tauStart[,g] <- sampleStats@means[[g]]
#     
#     # Current cov estimate:
#     curcov <- as.matrix(spectralshift(sampleStats@covs[[g]]))
#     
#     # Factor loadings:
#     for (f in seq_len(ncol(lambdaStart))){
#       # First principal component of sub cov:
#       if (any(lambda[,f,g]!=0)){
#         ev1 <- eigen(curcov[lambda[,f,g]!=0,lambda[,f,g]!=0])$vectors[,1]
#         lambdaStart[lambdaStart[,f,g]!=0,f,g] <- ev1 / ev1[1]        
#       } 
#     }
#     
#     # Residual variances, let's start by putting the vars on 1/4 times the observed variances:
#     Theta <- diag(diag(curcov)/4)
#     
#     # Check if this is positive definite:
#     ev <- eigen(curcov - Theta)$values
#     
#     # Shrink until it is positive definite:
#     loop <- 0
#     repeat{
#       ev <- eigen(curcov - Theta)$values
#       if (loop == 100){
#         # give up...
#  
#         Theta <- diag(nrow(Theta))
#         break
#       }
#       if (all(ev>0)){
#         break
#       }
#       Theta <- Theta * 0.9
#       loop <- loop + 1
#     }
# 
#     # Set delta to it's inversed square root:
#     # delta_epsilon_start[,,g] <- (1*delta_epsilon_start[,,g]!=0) * sqrt(corpcor::pseudoinverse(Theta))
# 
#     delta_epsilon_start[,,g] <- (1*delta_epsilon_start[,,g]!=0) * corpcor::pseudoinverse(sqrt(corpcor::pseudoinverse(Theta)))
#     
#     # This means that the factor-part is expected to be:
#     factorPart <- curcov - Theta
#     
#     # Let's take a pseudoinverse:
#     inv <- corpcor::pseudoinverse(kronecker(lambdaStart[,,g],lambdaStart[,,g]))
#     
#     # And obtain psi estimate:
#     PsiEstimate <- matrix(inv %*% as.vector(factorPart),nLatent,nLatent)
#     
#     # Let's put this as starting value:
#     sigma_zeta_start[,,g] <- (sigma_zeta[,,g]!=0) * PsiEstimate
#   }
# 
# 
#   # Generate the full parameter table:
#   pars <- generateAllParameterTables(
#     # Tau:
#     list(tau,
#          mat =  "tau",
#          op =  "~1",
#          symmetrical= FALSE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = "1",
#          start = tauStart),
#     
#     # mu_eta:
#     list(tau_eta,
#          mat =  "tau_eta",
#          op =  "~1",
#          symmetrical= FALSE, 
#          sampletable=sampleStats,
#          rownames = latents,
#          colnames = "1"),
#     
#     # Lambda:
#     list(lambda,
#          mat =  "lambda",
#          op =  "~=",
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = latents,
#          sparse = TRUE,
#          start = lambdaStart
#     ),
#     
#     # Beta:
#     list(beta,
#          mat =  "beta",
#          op =  "<-",
#          sampletable=sampleStats,
#          rownames = latents,
#          colnames = latents,
#          sparse = TRUE,
#          start = betaStart
#     ),
#     
#     # Latent covs:
#     list(sigma_zeta,
#          mat =  "sigma_zeta",
#          op =  "~~",
#          symmetrical= TRUE, 
#          sampletable=sampleStats,
#          rownames = latents,
#          colnames = latents,
#          sparse = TRUE,
#          posdef = TRUE,
#          start = sigma_zeta_start
#     ),
#     
#     # Reisdual network
#     list(omega_epsilon,
#          mat =  "omega_epsilon",
#          op =  "--",
#          symmetrical= TRUE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = sampleStats@variables$label,
#          sparse = TRUE,
#          posdef = TRUE,
#          diag0=TRUE
#     ),
#     
#     
#     # Delta_eta:
#     list(delta_epsilon,
#          mat =  "delta_epsilon",
#          op =  "~/~",
#          symmetrical= TRUE, 
#          sampletable=sampleStats,
#          rownames = sampleStats@variables$label,
#          colnames = sampleStats@variables$label,
#          sparse = TRUE,
#          posdef = TRUE,
#          diagonal = TRUE,
#          lower = 0,
#          start = delta_epsilon_start
#     )
#   )
# 
#   # Store in model:
#   model@parameters <- pars$partable
#   model@matrices <- pars$mattable
#   model@extramatrices <- list(
#     D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
#     Deta = psychonetrics::duplicationMatrix(nLatent), # non-strict duplciation matrix
#     L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
#     Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
#     In = Diagonal(nNode), # Identity of dim n
#     Inlatent = Diagonal(nLatent),
#     C = as(lavaan::lav_matrix_commutation(nNode, nLatent),"sparseMatrix"),
#     Cbeta = as(lavaan::lav_matrix_commutation(nLatent, nLatent),"sparseMatrix"),
#     A = psychonetrics::diagonalizationMatrix(nNode)
#   )
#     
#   
#   # Form the model matrices
#   model@modelmatrices <- formModelMatrices(model)
#   
#   
#   ### Baseline model ###
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
#                                                 lowertri = "full", 
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
#   # TODO: Identify model:
#   if (identify){
#     model <- identify(model)
#   }
#   # 
#   # Return model:
#   return(model)
# }