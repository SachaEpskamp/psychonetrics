# Latent network model creator
meta_varcov <- function(
  cors, # List of correlation matrices as input. Must contain NAs
  nobs, # vector of sample sizes as input
  
  Vmats, # Optional list of V matrices for each group. Will be averaged.
  # Vmethod = c("default","individual","weighted","psychonetrics_individual", "psychonetrics_weighted", "psychonetrics_pooled", "metaSEM_individual","metaSEM_weighted"), # How to obtain V matrices if Vmats is not supplied?
  # Vestimation = c("averaged","per_study"),

  Vmethod = c("individual","pooled","metaSEM_individual","metaSEM_weighted"), # How to obtain V matrices if Vmats is not supplied?
  Vestimation = c("averaged","per_study"),  
  
  # Model setup:
  type = c("cor", "ggm"), # Same as in varcov. Currently only cor and ggm are supported.
  sigma_y = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa_y = "full", # Precision
  # rho = "full", # Correlations
  omega_y = "full", # Partial correlations
  lowertri_y = "full", # Cholesky
  delta_y = "full", # Used for both ggm and pcor
  rho_y = "full", # Used for cor
  SD_y = "full", # Used for cor
  
  # Random effects setup:
  randomEffects = c("chol","cov","prec","ggm","cor"),
  sigma_randomEffects = "full", # (only lower tri is used) "empty", "full" or kappa structure, array (nvar * nvar * ngroup). NA indicates free, numeric indicates equality constraint, numeric indicates constraint
  kappa_randomEffects = "full", # Precision
  # rho = "full", # Correlations
  omega_randomEffects = "full", # Partial correlations
  lowertri_randomEffects = "full", # Cholesky
  delta_randomEffects = "full", # Used for both ggm and pcor
  rho_randomEffects = "full", # Used for cor
  SD_randomEffects = "full", # Used for cor
  
  vars, # Only used for labeling the variables
  
  # Some extra stuff:
  baseline_saturated = TRUE, # Leave to TRUE! Only used to stop recursive calls
  optimizer,
  estimator = c("FIML","ML"),
  
  sampleStats, # Leave to missing
  verbose = FALSE
){

  # warning("'meta_varcov' is still experimental.")
  sampleSizes <- nobs # FIXME
  
  # For now, I will always assume correlations were used in the input:
  corinput <- TRUE
  estimator <- match.arg(estimator)
  
  randomEffects <- match.arg(randomEffects)
  type <- match.arg(type)
  Vmethod <- match.arg(Vmethod)
  if (Vmethod == "default"){
    Vmethod <- "individual"
    # Vmethod <- 'psychonetrics_individual'
    if (!missing(Vmats)){
      stop("'Vmats' must be missing if 'Vmethod' is not 'default'")
    }
  }
  Vestimation <- match.arg(Vestimation)
  
  # set the labels:
  if (missing(vars)){
    vars <- unique(unlist(lapply(cors,colnames)))
    if (is.null(vars)){
      vars <- paste0("V",seq_len(max(sapply(cors,ncol))))
      if (length(unique(sapply(cors,colnames))) > 1){
        stop("Correlation matrices of different dimensions not supported without column labels.")
      }
      for (i in seq_along(cors)){
        rownames(cors[[i]]) <- colnames(cors[[i]]) <- vars
      }
    }
  } else {
    if (is.null(colnames(cors[[1]]))){
      for (i in seq_along(cors)){
        rownames(cors[[i]]) <- colnames(cors[[i]]) <- vars
      }
    }
  }
  
  # Number of nodes:
  nNode <- length(vars)
  
  # Reorder correlation matrices and add NAs:
  corsOld <- cors
  cors <- list()
  for (i in seq_along(corsOld)){
    cors[[i]] <- matrix(NA,nrow=nNode, ncol = nNode)
    rownames(cors[[i]]) <- colnames(cors[[i]]) <- vars
    varsOfStudy <- colnames(corsOld[[i]])
    matched <- match(varsOfStudy,vars)
    cors[[i]][matched,matched] <- as.matrix(corsOld[[i]])
  }
  
  # Form the dataset from the lower triangles of cor matrices:
  data <- dplyr::bind_rows(lapply(cors,function(x){
    df <- as.data.frame(t(x[lower.tri(x)]))
    names(df) <- paste0(rownames(x)[row(x)[lower.tri(x)]], " -- ",colnames(x)[col(x)[lower.tri(x)]])
    df
  }))
  
  # Labels for the correlations:
  corvars <- colnames(data)
  
  # Obtain sample stats:
  if (missing(sampleStats)){
    
    sampleStats <- samplestats(data = data,
                               vars = corvars,
                               missing = "pairwise",
                               fimldata = estimator == "FIML",
                               storedata = FALSE,
                               meanstructure = TRUE,
                               verbose=verbose,
                               fullFIML = (Vestimation == "per_study")
                               # fullFIML = (Vmethod == "individual")
                               )
  }
  
  # Overwrite corInput:
  sampleStats@corinput <- corinput
  
  
  # Check some things:
  nCor <- nrow(sampleStats@variables)
  
  # Generate model object:
  model <- generate_psychonetrics(model = "meta_varcov", sample = sampleStats, computed = FALSE,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  types = list(y = type, randomEffects = randomEffects),
                                  submodel = type, meanstructure = TRUE, verbose = verbose)
  
  # Number of groups:
  nGroup <- 1
  
  
  # Number of means and thresholds:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))
  
  # Add number of observations:
  model@sample@nobs <-
    nCor * (nCor-1) / 2 * nGroup + # Covariances per group
    nCor * nGroup + # Variances (ignored if correlation matrix is input)
    nMeans
  
  
  ### Estimate V matrices
  if (!missing(Vmats)){
    avgVmat <- Reduce("+", Vmats) / length(Vmats)
    
  } else {
    
    if (verbose){
      message("Computing sampling error approximation...")
    }
    
    
    # For the elimination matrix, I need this dummy matrix with indices:
    dumSig <- matrix(0,nNode,nNode)
    dumSig[lower.tri(dumSig,diag=FALSE)] <- seq_len(sum(lower.tri(dumSig,diag=FALSE)))
    
    # Every method should give Vmats and avgVmat...
    
    if (Vmethod == "individual"){ 
      # For each group, make a model and obtain VCOV:
      Vmats <- lapply(seq_along(cors),function(i){

        # Find the missing nodes:
        obs <- !apply(cors[[i]],2,function(x)all(is.na(x)))
        
        # Indices:
        inds <- c(dumSig[obs,obs,drop=FALSE])
        inds <- inds[inds!=0]
        
        # Elimintation matrix:
        L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),nNode*(nNode-1)/2))
        L <- as(L, "dMatrix")
        
        # Now obtain only the full subset correlation matrix:
        cmat <- as(cors[[i]][obs,obs], "matrix")
        
        k <- solve_symmetric_cpp_matrixonly(cmat)
        D2 <- duplicationMatrix(ncol(cmat), FALSE)
        v <- 0.5 * nobs[i] * t(D2) %*% (k %x% k) %*% D2
        vcov <- solve_symmetric_cpp_matrixonly(as.matrix(v))
        
        
        
        # Now expand using the elmination matrix:
        res <- as.matrix(t(L) %*% vcov %*% L)
        return(0.5 * (res + t(res)))
      })
      
      avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))
      
    } else {
      
        # If there are any NAs, use covariance input to compute model using FIML:
        if(any(is.na(unlist(lapply(cors,as.vector))))){
          # browser()
          # Single multi-group model:
          mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = c("SD","rho","mu"), baseline_saturated = FALSE, verbose = FALSE,
                        estimator = "FIML", covtype = "ML", meanstructure = TRUE)
          mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
          acov <- getVCOV(mod)
          avgVmat <- acov
          ind <- which(mod@parameters$matrix[match(seq_len(max(mod@parameters$par)),mod@parameters$par)] == "rho")
          avgVmat <- acov[ind,ind] * length(cors)

        } else {
          mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = "rho", baseline_saturated = FALSE, verbose = FALSE)
          mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
          acov <- getVCOV(mod)
          avgVmat <- acov * length(cors)
        }

        # Compute Vmat per dataset:
        Vmats <- lapply(nobs,function(n) mean(nobs)/n * avgVmat)
    }
    
    # 
    # if (Vmethod %in% c("weighted","psychonetrics_weighted")){
    #   
    #   # Make an average correlation matrix, weighted per sample size:
    #   avgCormat <- diag(1, nNode)
    #   for (i in 2:nNode){
    #     for (j in 1:(i-1)){
    #       avgCormat[i,j] <- avgCormat[j,i] <- weighted.mean(sapply(cors,'[',i,j), nobs, na.rm = TRUE)
    #     }
    #   }
    #   avgN <- mean(nobs,na.rm=TRUE)
    # }
    # 
    # 
    # if (Vmethod == "weighted"){
    #   
    #   k <- solve(avgCormat)
    #   D2 <- duplicationMatrix(ncol(avgCormat), FALSE)
    #   v <- 0.5 * avgN * t(D2) %*% (k %x% k) %*% D2
    #   avgVmat <- solve(v)
    #   
    #   # Compute Vmat per dataset:
    #   Vmats <- lapply(nobs,function(n) mean(nobs)/n * avgVmat)
    # }
    # 
    # if (Vmethod == "psychonetrics_weighted"){
    #   
    #   mod <- varcov(covs=avgCormat,nobs=avgN, corinput = TRUE, type = "cor", baseline_saturated = FALSE, verbose = FALSE)
    #   mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
    #   avgVmat <- getVCOV(mod)
    #   
    #   
    #   # Compute Vmat per dataset:
    #   Vmats <- lapply(nobs,function(n) mean(nobs)/n * avgVmat)
    # }
    # 
    # if (Vmethod == "psychonetrics_individual"){
    #   
    #   # For each group, make a model and obtain VCOV:
    #   Vmats <- lapply(seq_along(cors),function(i){
    #     
    #     # Find the missing nodes:
    #     obs <- !apply(cors[[i]],2,function(x)all(is.na(x)))
    #     
    #     # Indices:
    #     inds <- c(dumSig[obs,obs,drop=FALSE])
    #     inds <- inds[inds!=0]
    #     
    #     # Elimintation matrix:
    #     L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),nNode*(nNode-1)/2))
    #     L <- as(L, "indMatrix")
    #     
    #     # Now obtain only the full subset correlation matrix:
    #     cmat <- as(cors[[i]][obs,obs], "matrix")
    #     
    #     # Now run psychonetrics:
    #     mod <- varcov(covs=cmat,nobs=sampleSizes[i], corinput = TRUE, type = "cor", baseline_saturated = FALSE, verbose = FALSE)
    #     mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
    #     vcov <- getVCOV(mod)
    #     
    #     # Now expand using the elmination matrix:
    #     t(L) %*% vcov %*% L
    #     
    #     
    #     
    #   })
    #   
    #   avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))
    #   
    # }
    # 
    # 
    # if (Vmethod == "psychonetrics_pooled"){
    #   
    #   
    #   # If there are any NAs, use covariance input to compute model using FIML:
    #   if(any(is.na(unlist(lapply(cors,as.vector))))){
    #     # Single multi-group model:
    #     mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = c("SD","rho","mu"), baseline_saturated = FALSE, verbose = FALSE,
    #                   estimator = "FIML", covtype = "ML", meanstructure = TRUE)
    #     mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
    #     acov <- getVCOV(mod)
    #     
    #     ind <- which(mod@parameters$matrix[match(seq_len(max(mod@parameters$par)),mod@parameters$par)] == "rho")
    #     avgVmat <- acov[ind,ind] * length(cors)
    #     
    #   } else {
    #     mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = "rho", baseline_saturated = FALSE, verbose = FALSE)        
    #     mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
    #     acov <- getVCOV(mod)
    #     avgVmat <- acov * length(cors)
    #   }
    #   
    #   # Compute Vmat per dataset:
    #   Vmats <- lapply(nobs,function(n) mean(nobs)/n * avgVmat)
    #   
    #   # 
    #   # # mod <- varcov(covs=cors,nobs=sampleSizes, corinput = TRUE, type =  "cor", equal = "rho", baseline_saturated = FALSE, verbose = FALSE,
    #   #               # estimator = "FIML")
    #   # mod <- runmodel(mod, addfit = FALSE, addMIs = FALSE, addSEs = FALSE, verbose = FALSE)
    #   # acov <- getVCOV(mod)
    #   # avgVmat <- acov * length(cors)
    # }
    # 
    if (Vmethod == "metaSEM_individual"){

      acovs <- metaSEM::asyCov(cors, sampleSizes, acov = "individual")
      acovs[is.na(acovs)] <- 0
      Vmats <- list()
      for (i in seq_len(nrow(acovs))){
        Vmats[[i]] <- matrix(0, nCor, nCor)
        Vmats[[i]][lower.tri(Vmats[[i]],diag=TRUE)] <- acovs[i,]
        Vmats[[i]][upper.tri(Vmats[[i]],diag=TRUE)] <- t(Vmats[[i]])[upper.tri(Vmats[[i]],diag=TRUE)]
      }
      avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))

    }



    if (Vmethod == "metaSEM_weighted"){
      acovs <- metaSEM::asyCov(cors, sampleSizes, acov = "weighted")
      acovs[is.na(acovs)] <- 0
      Vmats <- list()
      for (i in seq_len(nrow(acovs))){
        Vmats[[i]] <- matrix(0, nCor, nCor)
        Vmats[[i]][lower.tri(Vmats[[i]],diag=TRUE)] <- acovs[i,]
        Vmats[[i]][upper.tri(Vmats[[i]],diag=TRUE)] <- t(Vmats[[i]])[upper.tri(Vmats[[i]],diag=TRUE)]
      }
      avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))
    }
  }
  
  ####
  
  
  # Model matrices:
  modMatrices <- list()
  
  # Obtain the expected homogeneous cor structure from means:
  expCorsVec <- model@sample@means[[1]]
  expCors <- matrix(1,nNode,nNode)
  expCors[lower.tri(expCors)] <- expCorsVec
  expCors[upper.tri(expCors)] <- t(expCors)[upper.tri(expCors)]
  
  
  # Fix sigma
  if (type == "cov"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'cov'. Use type = 'cor' or set corinput = FALSE")
    }
    
    modMatrices$sigma_y <- matrixsetup_sigma(sigma_y,
                                             expcov=list(expCors),
                                             nNode = nNode,
                                             nGroup = nGroup,
                                             labels = vars,
                                             equal = FALSE,
                                             sampletable = sampleStats,
                                             name = "sigma_y")
  } else if (type == "chol"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'chol'.")
    }
    modMatrices$lowertri_y <- matrixsetup_lowertri(lowertri_y,
                                                   expcov=list(expCors),
                                                   nNode = nNode,
                                                   nGroup = nGroup,
                                                   labels = vars,
                                                   equal = FALSE,
                                                   sampletable = sampleStats,
                                                   name = "lowertri_y")
  } else if (type == "ggm"){
    # Add omega matrix:
    modMatrices$omega_y <- matrixsetup_omega(omega_y,
                                             expcov=list(expCors),
                                             nNode = nNode,
                                             nGroup = nGroup,
                                             labels = vars,
                                             equal = FALSE,
                                             sampletable = sampleStats,
                                             name = "omega_y")
    
    if (!corinput){
      # Add delta matrix (ingored if corinput == TRUE):
      modMatrices$delta_y <- matrixsetup_delta(delta_y,
                                               expcov=list(expCors),
                                               nNode = nNode,
                                               nGroup = nGroup,
                                               labels = vars,
                                               equal = FALSE,
                                               sampletable = sampleStats,
                                               name = "delta_y",
                                               omegaStart =  modMatrices$omega_y$start)
    }
    
  } else if (type == "prec"){
    if (corinput){
      stop("Correlation matrix input is not supported for type = 'prec'. Use type = 'ggm' or set corinput = FALSE")
    }
    
    # Add omega matrix:
    modMatrices$kappa_y <- matrixsetup_kappa(kappa_y,
                                             expcov=list(expCors),
                                             nNode = nNode,
                                             nGroup = nGroup,
                                             labels = vars,
                                             equal = FALSE,
                                             sampletable = sampleStats,
                                             name = "kappa_y")
  } else if (type == "cor"){
    # Add rho matrix:
    modMatrices$rho_y <- matrixsetup_rho(rho_y,
                                         expcov=list(expCors),
                                         nNode = nNode,
                                         nGroup = nGroup,
                                         labels = vars,
                                         equal = FALSE,
                                         sampletable = sampleStats,
                                         name = "rho_y")
    
    if (!corinput){
      # Add SD matrix (ignored if corinput == TRUE):
      modMatrices$SD_y <- matrixsetup_SD(SD_y,
                                         expcov=list(expCors),
                                         nNode = nNode,
                                         nGroup = nGroup,
                                         labels = vars,
                                         equal = FALSE,
                                         sampletable = sampleStats,
                                         name = "rho_y")
    }
  }
  
  
  # Compute expected random effects matrix:
  expRanEffects <- as.matrix(spectralshift(sampleStats@covs[[1]] - avgVmat))
  
  
  # Add random effets matrices:
  if (randomEffects == "cov"){
    modMatrices$sigma_randomEffects <- matrixsetup_sigma(sigma_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "sigma_randomEffects")
  } else if (randomEffects == "chol"){
    
    modMatrices$lowertri_randomEffects <- matrixsetup_lowertri(lowertri_randomEffects,
                                                               expcov=list(expRanEffects),
                                                               nNode = nCor,
                                                               nGroup = nGroup,
                                                               labels = corvars,
                                                               equal = FALSE,
                                                               sampletable = sampleStats,
                                                               name = "lowertri_randomEffects")
  } else if (randomEffects == "ggm"){
    # Add omega matrix:
    modMatrices$omega_randomEffects <- matrixsetup_omega(omega_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "omega_randomEffects")
    
    
    # Add delta matrix (ingored if corinput == TRUE):
    modMatrices$delta_randomEffects <- matrixsetup_delta(delta_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "delta_randomEffects",
                                                         omegaStart =  modMatrices$omega_randomEffects$start)
    
    
  } else if (randomEffects == "prec"){
    
    # Add omega matrix:
    modMatrices$kappa_randomEffects <- matrixsetup_kappa(kappa_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "kappa_randomEffects")
  } else if (randomEffects == "cor"){
    # Add rho matrix:
    modMatrices$rho_randomEffects <- matrixsetup_rho(rho_randomEffects,
                                                     expcov=list(expRanEffects),
                                                     nNode = nCor,
                                                     nGroup = nGroup,
                                                     labels = corvars,
                                                     equal = FALSE,
                                                     sampletable = sampleStats,
                                                     name = "rho_randomEffects")
    
    
    # Add SD matrix (ignored if corinput == TRUE):
    modMatrices$SD_randomEffects <- matrixsetup_SD(SD_randomEffects,
                                                   expcov=list(expRanEffects),
                                                   nNode = nCor,
                                                   nGroup = nGroup,
                                                   labels = corvars,
                                                   equal = FALSE,
                                                   sampletable = sampleStats,
                                                   name = "SD_randomEffects")
  }
  
  
  
  # Generate the full parameter table:
  pars <- do.call(generateAllParameterTables, modMatrices)
  
  
  # Store in model:
  model@parameters <- pars$partable
  model@matrices <- pars$mattable
  
  
  # Just add all the matrices
  # FIXME: This stores some objects that may not be needed.
  model@extramatrices <- list(
    D = psychonetrics::duplicationMatrix(nNode), # non-strict duplciation matrix
    L = psychonetrics::eliminationMatrix(nNode), # Elinimation matrix
    Lstar = psychonetrics::eliminationMatrix(nNode, diag=FALSE), # Elinimation matrix
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE), # Strict duplicaton matrix
    In = as(diag(nNode),"dMatrix"), # Identity of dim n
    A = psychonetrics::diagonalizationMatrix(nNode),
    C = as(lavaan::lav_matrix_commutation(nNode,nNode),"dMatrix"),
    
    # Random effects:
    D_c = psychonetrics::duplicationMatrix(nCor), # non-strict duplciation matrix
    L_c = psychonetrics::eliminationMatrix(nCor), # Elinimation matrix
    Lstar_c = psychonetrics::eliminationMatrix(nCor, diag=FALSE), # Elinimation matrix
    Dstar_c = psychonetrics::duplicationMatrix(nCor,diag = FALSE), # Strict duplicaton matrix
    In_c = as(diag(nCor),"dMatrix"), # Identity of dim n
    A_c = psychonetrics::diagonalizationMatrix(nCor),
    C_c = as(lavaan::lav_matrix_commutation(nCor,nCor),"dMatrix"),
    
    # Add the vmat:
    V = avgVmat,
    Vall = Vmats,
    Vmethod = Vmethod,
    Vestimation = Vestimation # FIXME: This is nicer somewhere else...
  )
  

  
  # 
  # if (Vmethod == "pooled"){
  #   model@extramatrices$Vall <- avgVmat
  # } else {
  #   model@extramatrices$Vmats <- Vmats
  # }
  # 
  
  
  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)
  
  ### Baseline model ###
  if (baseline_saturated){
    
    # FIXME: Dummy sample stats:
    sampleStats2 <- sampleStats
    sampleStats2@corinput <- FALSE
    
    # Form baseline model:
    model@baseline_saturated$baseline <- varcov(data,
                                                mu = rep(0,nCor),
                                                type = "chol",
                                                lowertri = "empty",
                                                vars = corvars,
                                                missing = missing,
                                                estimator = estimator,
                                                baseline_saturated = FALSE,
                                                sampleStats=sampleStats2)
    
    
    
    model@baseline_saturated$baseline@sample@fullFIML <- FALSE
    # Add model:
    # model@baseline_saturated$baseline@fitfunctions$extramatrices$M <- Mmatrix(model@baseline_saturated$baseline@parameters)
    
    
    ### Saturated model ###
    
    model@baseline_saturated$saturated <- varcov(data,
                                                 type = "chol",
                                                 lowertri = "full",
                                                 vars = corvars,
                                                 missing = missing,
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats=sampleStats2)
    
    model@baseline_saturated$saturated@sample@fullFIML <- FALSE
    
    # if not FIML, Treat as computed:
    if (estimator != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      
      # FIXME: TODO
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
    }
  }
  
  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }
  
  # Return model:
  return(model)
}
