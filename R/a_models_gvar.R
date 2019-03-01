# Latent network model creator
gvar <- function(
  data, # Dataset
  beta = "full",
  omega_zeta = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  delta_zeta = "full", # If missing, just full for both groups or equal
  mu,
  beepvar,
  dayvar,
  idvar,
  vars, # character indicating the variables Extracted if missing from data - group variable
  groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
  covs, # alternative covs (array nvar * nvar * ngroup)
  means, # alternative means (matrix nvar * ngroup)
  nobs, # Alternative if data is missing (length ngroup)
  missing = "fiml",
  equal = "none", # Can also be any of the matrices
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  # fitfunctions, # Leave empty
  estimator = "ML",
  optimizer = "ucminf"
){
  # FIXME: Not sure why needed...
  if (missing(vars)) vars2 <- NULL else vars2 <- vars
  if (missing(idvar)) idvar <- NULL
  if (missing(dayvar)) dayvar <- NULL
  if (missing(beepvar)) beepvar <- NULL
  if (missing(groups)) groups <- NULL
  
  # If data is not missing, make augmented data:
  if (!missing(data)){
    data <- tsData(data, vars = vars2, beepvar = beepvar, dayvar = dayvar, idvar = idvar, groupvar = groups)
    if (is.null(groups)){
      vars <- colnames(data)  
    } else {
      vars <- colnames(data)[colnames(data)!=groups]
    }
    
  }
  # Obtain sample stats:
  sampleStats <- samplestats(data = data, 
                             vars = vars, 
                             groups = groups,
                             covs = covs, 
                             means = means, 
                             nobs = nobs, 
                             missing = missing)
  # Check if number of variables is an even number:
  if (nrow(sampleStats@variables)%%2!=0){
    stop("Number of variables is not an even number: variance-covariance matrix cannot be a Toeplitz matrix. ")
  }
  
  # Check some things:
  nNode <- nrow(sampleStats@variables) / 2
  
  # Generate model object:
  model <- generate_psychonetrics(model = "gvar",sample = sampleStats,computed = FALSE, 
                                  equal = equal,
                                  optimizer = optimizer, estimator = estimator, distribution = "Gaussian")
  
  # Number of groups:
  nGroup <- nrow(model@sample@groups)
  
  # Add number of observations:
  model@sample@nobs <-  
    nNode * (nNode+1) / 2 * nGroup + # Covariances per group
    nNode * nGroup # Means per group
  
  # Fix mu
  mu <- fixMu(mu,nGroup,nNode*2,"mu" %in% equal)
  
  # Fix the omega_zeta matrix:
  omega_zeta <- fixAdj(omega_zeta,nGroup,nNode,"omega_zeta" %in% equal,diag0=TRUE)
  
  # Fix delta:
  delta_zeta <- fixAdj(delta_zeta,nGroup,nNode,"delta_zeta" %in% equal,diagonal=TRUE)
  
  # Fix beta matrix:
  beta <- fixMatrix(beta,nGroup,nNode,nNode,"beta" %in% equal)
 
  # Exogenous varcov:
  exoVarCov <- fixAdj("full",nGroup,nNode)
  
  # Generate starting values:
  omegaStart <- omega_zeta
  deltaStart <- delta_zeta
  muStart <- mu
  betaStart <- beta
  exoStart <- exoVarCov
  for (g in 1:nGroup){
    # Means with sample means:
    muStart[,g] <- sampleStats@means[[g]]
    
    # Let's get three blocks:
    exoStart[,,g] <- as.matrix(sampleStats@covs[[g]][1:nNode,1:nNode])
    Sstar <- exoStart[,,g]
    S1 <- sampleStats@covs[[g]][nNode + (1:nNode),1:nNode]
    S0 <- sampleStats@covs[[g]][nNode + (1:nNode),nNode + (1:nNode)]
    S0inv <- corpcor::pseudoinverse(S0)
    
    # A prior guess for beta is:
    betaStart[,,g] <- (beta[,,g]!=0) * as.matrix(S1 %*% S0inv)
    
    # A prior guess for the contemporaneous covariances is (Schur complement):
    contCov <- as.matrix(Sstar - t(S1) %*% S0inv %*% S1)
    
    # Let's use this as starting estimate:
    zeroes <- which(omega_zeta[,,g]==0 & t(omega_zeta[,,g])==0 & diag(nNode) != 1,arr.ind=TRUE)
    if (nrow(zeroes) == 0){
      # glas <- glasso(as.matrix(t(lambdaStart[,,g]) %*% sampleStats@covs[[g]] %*% lambdaStart[,,g]),
      #                rho = 1e-10, zero = zeroes)      
      # wi <- corpcor::pseudoinverse(contCov)
      wi <- glasso(contCov, rho = 0.1)$wi
    } else {
      glas <- glasso(contCov,
                     rho = 0.1, zero = zeroes)
      wi <- glas$wi
    }
    
    # Network starting values:
    pcors <- as.matrix(qgraph::wi2net(wi))
    pcors[upper.tri(pcors)] <- 0
    omegaStart[,,g] <- ifelse(pcors!=0,pcors,ifelse(omegaStart[,,g]!=0,0.1,0))
    # omegaStart[,,g] <- 1*(omegaStart[,,g]!=0) * 0.05
    diag(omegaStart[,,g] ) <- 0
    # Delta:
    # deltaStart[,,g] <- diag(1/sqrt(diag(wi)))
    deltaStart[,,g] <- diag(sqrt(diag(wi)))
  }

  # Generate the full parameter table:
  pars <- generateAllParameterTables(
    # mu:
    list(mu,
         mat =  "mu",
         op =  "~1",
         symmetrical= FALSE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label,
         colnames = "1",
         start = muStart),
    
    # Exogenous variancecovariances:
    list(exoVarCov,
         mat =  "exogenous_sigma",
         op =  "~~",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label[1:nNode],
         colnames = sampleStats@variables$label[1:nNode],
         start = exoStart),
    
    # Beta:
    list(beta,
         mat =  "beta",
         op =  "<-",
         sampletable=sampleStats,
         rownames = sampleStats@variables$label[nNode + (1:nNode)],
         colnames = sampleStats@variables$label[1:nNode],
         sparse = TRUE,
         start = betaStart
    ),
    
    # Omega:
    list(omega_zeta,
         mat =  "omega_zeta",
         op =  "--",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label[nNode + (1:nNode)],
         colnames = sampleStats@variables$label[nNode + (1:nNode)],
         sparse = TRUE,
         posdef = TRUE,
         diag0=TRUE,
         lower = -1,
         upper = 1,
         start = omegaStart
    ),
    
    # Delta_zeta:
    list(delta_zeta,
         mat =  "delta_zeta",
         op =  "~/~",
         symmetrical= TRUE, 
         sampletable=sampleStats,
         rownames = sampleStats@variables$label[nNode + (1:nNode)],
         colnames = sampleStats@variables$label[nNode + (1:nNode)],
         sparse = TRUE,
         posdef = TRUE,
         diagonal = TRUE,
         lower = 0,
         start = deltaStart
    )
    
  )
  # 
  # # Set lambda start:
  # pars$partable$est[pars$partable$matrix=="lambda" & !pars$partable$fixed] <- 0.5
  # 
  # # Set mu start
  # for (g in 1:nGroup){
  #   pars$partable$est[pars$partable$matrix=="mu" & pars$partable$group_id == g] <- sampleStats@means[[g]]
  # }
  # 
  # Come up with P...
  # Dummy matrix to contain indices:
  dummySigma <- matrix(0,nNode*2,nNode*2)
  smallMat <- matrix(0,nNode,nNode)
  dummySigma[1:nNode,1:nNode][lower.tri(smallMat,diag=TRUE)] <- seq_len(nNode*(nNode+1)/2)
  dummySigma[nNode + (1:nNode),nNode + (1:nNode)][lower.tri(smallMat,diag=TRUE)] <- max(dummySigma) + seq_len(nNode*(nNode+1)/2)
  dummySigma[nNode + (1:nNode),1:nNode] <- max(dummySigma) + seq_len(nNode^2)
  inds <- dummySigma[lower.tri(dummySigma,diag=TRUE)]
  
  # P matrix:
  # P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=inds))
  P <- bdiag(Diagonal(nNode*2),sparseMatrix(j=seq_along(inds),i=order(inds)))

  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  model@extramatrices <- list(
    D =  psychonetrics::duplicationMatrix(nNode*2), # Toeplitz matrix D 
    D2 = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    In = Diagonal(nNode), # Identity of dim n
    In2 = Diagonal(nNode), # Identity of dim n^2
    A = psychonetrics::diagonalizationMatrix(nNode),
    P=P # Permutation matrix
  )
    
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  
  ### Baseline model ###
  if (baseline_saturated){
    # Via ggm:
    warning("Baseline and saturated are not correct yet")
    model@baseline_saturated$baseline <- ggm(data, 
                                             omega = "empty",
                                             vars = vars,
                                             groups = groups,
                                             covs = covs,
                                             means = means,
                                             nobs = nobs,
                                             missing = missing,
                                             equal = equal,
                                             baseline_saturated = FALSE)
    
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    model@baseline_saturated$saturated <- ggm(data, 
                                              omega = "full", 
                                              vars = vars,
                                              groups = groups,
                                              covs = covs,
                                              means = means,
                                              nobs = nobs,
                                              missing = missing,
                                              equal = equal,
                                              baseline_saturated = FALSE)
    
    # Treat as computed:
    model@baseline_saturated$saturated@computed <- TRUE
    
    # Add saturated fit
    model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
  }
  
 
  # Return model:
  return(model)
}
