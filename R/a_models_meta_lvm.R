# Meta-analytic latent variable model (single-stage random effects)
meta_lvm <- function(
  cors, # List of correlation matrices as input. Must contain NAs for missing variables
  nobs, # vector of sample sizes as input

  Vmats, # Optional list of V matrices for each group. Will be averaged.
  Vmethod = c("individual","pooled","metaSEM_individual","metaSEM_weighted"), # How to obtain V matrices if Vmats is not supplied?
  Vestimation = c("averaged","per_study"),

  # LVM structure:
  lambda, # REQUIRED: factor loading matrix (nvar x nlat), or a pattern matrix
  beta = "zero", # latent regression matrix

  # Latent covariance structure:
  latent = c("cov","chol","prec","ggm"),
  sigma_zeta = "full",
  kappa_zeta = "full",
  omega_zeta = "full",
  lowertri_zeta = "full",
  delta_zeta = "full",

  # Residual covariance structure:
  residual = c("cov","chol","prec","ggm"),
  sigma_epsilon = "diag",
  kappa_epsilon = "diag",
  omega_epsilon = "zero",
  lowertri_epsilon = "diag",
  delta_epsilon = "diag",

  # Identification:
  identification = "variance", # Variance identification for meta-analytic models

  # Random effects setup:
  randomEffects = c("chol","cov","prec","ggm","cor"),
  sigma_randomEffects = "full",
  kappa_randomEffects = "full",
  omega_randomEffects = "full",
  lowertri_randomEffects = "full",
  delta_randomEffects = "full",
  rho_randomEffects = "full",
  SD_randomEffects = "full",

  # Naming:
  vars, # character vector of variable names
  latents, # Name of latent variables

  # Some extra stuff:
  baseline_saturated = TRUE,
  optimizer,
  estimator = c("FIML","ML"),

  sampleStats,
  verbose = FALSE,
  bootstrap = FALSE,
  boot_sub,
  boot_resample
){

  message(paste0("Note: 'meta_lvm()' is experimental in psychonetrics ",
                 utils::packageVersion("psychonetrics"),
                 ". Please report any unexpected behavior to https://github.com/SachaEpskamp/psychonetrics/issues"))

  sampleSizes <- nobs

  # Treat input correlations as covariances (with free sigma_epsilon diagonal):
  corinput <- FALSE
  estimator <- match.arg(estimator)

  randomEffects <- match.arg(randomEffects)
  latent <- match.arg(latent)
  residual <- match.arg(residual)
  Vmethod <- match.arg(Vmethod)
  Vestimation <- match.arg(Vestimation)

  # lambda is required:
  if (missing(lambda)){
    stop("'lambda' may not be missing")
  }
  if (is.character(lambda)){
    stop("'lambda' may not be a string")
  }

  # Set the labels:
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
                               fullFIML = (Vestimation == "per_study"),
                               bootstrap=bootstrap,
                               boot_sub = boot_sub,
                               boot_resample = boot_resample)
  }

  # Treat correlations as covariances (no corinput constraint):
  sampleStats@corinput <- FALSE

  # Check some things:
  nCor <- nrow(sampleStats@variables)

  # Number of latents:
  nLatent <- ncol(lambda)

  # If latents is not provided, make it:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }

  # Generate model object:
  model <- generate_psychonetrics(model = "meta_lvm", sample = sampleStats, computed = FALSE,
                                  optimizer =  defaultoptimizer(), estimator = estimator, distribution = "Gaussian",
                                  types = list(latent = latent, residual = residual, randomEffects = randomEffects),
                                  meanstructure = TRUE, verbose = verbose)

  # Number of groups:
  nGroup <- 1

  # Number of means:
  nMeans <- sum(sapply(model@sample@means,function(x)sum(!is.na(x))))

  # Add number of observations:
  model@sample@nobs <-
    nCor * (nCor-1) / 2 * nGroup + # Covariances per group
    nCor * nGroup + # Variances
    nMeans

  ### Estimate V matrices ###
  if (!missing(Vmats)){
    avgVmat <- Reduce("+", Vmats) / length(Vmats)
  } else {
    if (verbose){
      message("Computing sampling error approximation...")
    }

    # For the elimination matrix, I need this dummy matrix with indices:
    dumSig <- matrix(0,nNode,nNode)
    dumSig[lower.tri(dumSig,diag=FALSE)] <- seq_len(sum(lower.tri(dumSig,diag=FALSE)))

    if (Vmethod == "individual"){
      # For each group, make a model and obtain VCOV:
      Vmats <- lapply(seq_along(cors),function(i){
        # Find the missing nodes:
        obs <- !apply(cors[[i]],2,function(x)all(is.na(x)))

        # Indices:
        inds <- c(dumSig[obs,obs,drop=FALSE])
        inds <- inds[inds!=0]

        # Elimination matrix:
        L <- sparseMatrix(i=seq_along(inds),j=inds,dims=c(length(inds),nNode*(nNode-1)/2))
        L <- as(L, "dMatrix")

        # Now obtain only the full subset correlation matrix:
        cmat <- as(cors[[i]][obs,obs], "matrix")

        k <- solve_symmetric_cpp_matrixonly(cmat)
        D2 <- duplicationMatrix(ncol(cmat), FALSE)
        v <- 0.5 * nobs[i] * t(D2) %*% (k %x% k) %*% D2
        vcov <- solve_symmetric_cpp_matrixonly(as.matrix(v))

        # Now expand using the elimination matrix:
        res <- as.matrix(t(L) %*% vcov %*% L)
        return(0.5 * (res + t(res)))
      })

      avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))

    } else if (Vmethod == "pooled") {
      # If there are any NAs, use covariance input to compute model using FIML:
      if(any(is.na(unlist(lapply(cors,as.vector))))){
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

    } else if (Vmethod == "metaSEM_individual"){
      acovs <- metaSEM::asyCov(cors, sampleSizes, acov = "individual")
      acovs[is.na(acovs)] <- 0
      Vmats <- list()
      for (i in seq_len(nrow(acovs))){
        Vmats[[i]] <- matrix(0, nCor, nCor)
        Vmats[[i]][lower.tri(Vmats[[i]],diag=TRUE)] <- acovs[i,]
        Vmats[[i]][upper.tri(Vmats[[i]],diag=TRUE)] <- t(Vmats[[i]])[upper.tri(Vmats[[i]],diag=TRUE)]
      }
      avgVmat <- Reduce("+", Vmats) / Reduce("+",lapply(Vmats,function(x)x!=0))

    } else if (Vmethod == "metaSEM_weighted"){
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

  #### LVM Model matrices ####
  modMatrices <- list()

  # Setup lambda (using dummy covariance matrices for starting values):
  # Construct dummy covariance matrices from the average correlation:
  expCorsVec <- model@sample@means[[1]]
  expCors <- matrix(1,nNode,nNode)
  expCors[lower.tri(expCors)] <- expCorsVec
  expCors[upper.tri(expCors)] <- t(expCors)[upper.tri(expCors)]

  modMatrices$lambda <- matrixsetup_lambda(lambda, expcov=list(expCors), nGroup = nGroup, equal = FALSE,
                                           observednames = vars, latentnames = latents,
                                           sampletable = sampleStats, identification = identification, simple = FALSE)

  # Compute the expected latent and residual cov matrices from lambda setup:
  expLatSigma <- lapply(1:nGroup,function(x)matrix(0,nLatent,nLatent))
  expResidSigma <- lapply(1:nGroup,function(x)matrix(0,nNode,nNode))

  for (g in 1:nGroup){
    expResidSigma[[g]] <- modMatrices$lambda$sigma_epsilon_start[,,g]
    expLatSigma[[g]] <-  modMatrices$lambda$sigma_zeta_start[,,g]
  }

  # Setup beta:
  modMatrices$beta <- matrixsetup_beta(beta, nNode = nLatent, nGroup = nGroup,
                                        labels = latents, sampletable = sampleStats,
                                        equal = FALSE)

  # Latent varcov:
  if (latent == "cov"){
    modMatrices$sigma_zeta <- matrixsetup_sigma(sigma_zeta,
                                                name = "sigma_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = FALSE, sampletable = sampleStats,
                                                beta = modMatrices$beta[[1]])
  } else if (latent == "chol"){
    modMatrices$lowertri_zeta <- matrixsetup_lowertri(lowertri_zeta,
                                                      name = "lowertri_zeta",
                                                      expcov=expLatSigma,
                                                      nNode = nLatent,
                                                      nGroup = nGroup,
                                                      labels = latents,
                                                      equal = FALSE, sampletable = sampleStats,
                                                      beta = modMatrices$beta[[1]])
  } else if (latent == "ggm"){
    modMatrices$omega_zeta <- matrixsetup_omega(omega_zeta,
                                                name = "omega_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = FALSE, sampletable = sampleStats,
                                                beta = modMatrices$beta[[1]])
    modMatrices$delta_zeta <- matrixsetup_delta(delta_zeta,
                                                name = "delta_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = FALSE, sampletable = sampleStats,
                                                omegaStart =  modMatrices$omega_zeta$start)
  } else if (latent == "prec"){
    modMatrices$kappa_zeta <- matrixsetup_kappa(kappa_zeta,
                                                name = "kappa_zeta",
                                                expcov=expLatSigma,
                                                nNode = nLatent,
                                                nGroup = nGroup,
                                                labels = latents,
                                                equal = FALSE, sampletable = sampleStats,
                                                beta = modMatrices$beta[[1]])
  }

  ### Residual varcov ###
  if (residual == "cov"){
    modMatrices$sigma_epsilon <- matrixsetup_sigma(sigma_epsilon,
                                                   name = "sigma_epsilon",
                                                   expcov=expResidSigma,
                                                   nNode = nNode,
                                                   nGroup = nGroup,
                                                   labels = vars,
                                                   equal = FALSE, sampletable = sampleStats)
  } else if (residual == "chol"){
    modMatrices$lowertri_epsilon <- matrixsetup_lowertri(lowertri_epsilon,
                                                         name = "lowertri_epsilon",
                                                         expcov=expResidSigma,
                                                         nNode = nNode,
                                                         nGroup = nGroup,
                                                         labels = vars,
                                                         equal = FALSE, sampletable = sampleStats)
  } else if (residual == "ggm"){
    modMatrices$omega_epsilon <- matrixsetup_omega(omega_epsilon,
                                                   name = "omega_epsilon",
                                                   expcov=expResidSigma,
                                                   nNode = nNode,
                                                   nGroup = nGroup,
                                                   labels = vars,
                                                   equal = FALSE, sampletable = sampleStats)
    modMatrices$delta_epsilon <- matrixsetup_delta(delta_epsilon,
                                                   name = "delta_epsilon",
                                                   expcov=expResidSigma,
                                                   nNode = nNode,
                                                   nGroup = nGroup,
                                                   labels = vars,
                                                   equal = FALSE, sampletable = sampleStats,
                                                   omegaStart =  modMatrices$omega_epsilon$start)
  } else if (residual == "prec"){
    modMatrices$kappa_epsilon <- matrixsetup_kappa(kappa_epsilon,
                                                   name = "kappa_epsilon",
                                                   expcov=expResidSigma,
                                                   nNode = nNode,
                                                   nGroup = nGroup,
                                                   labels = vars,
                                                   equal = FALSE, sampletable = sampleStats)
  }

  #### Random effects matrices ####
  # Compute expected random effects matrix:
  expRanEffects <- as.matrix(spectralshift(sampleStats@covs[[1]] - avgVmat))

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
    modMatrices$omega_randomEffects <- matrixsetup_omega(omega_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "omega_randomEffects")
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
    modMatrices$kappa_randomEffects <- matrixsetup_kappa(kappa_randomEffects,
                                                         expcov=list(expRanEffects),
                                                         nNode = nCor,
                                                         nGroup = nGroup,
                                                         labels = corvars,
                                                         equal = FALSE,
                                                         sampletable = sampleStats,
                                                         name = "kappa_randomEffects")
  } else if (randomEffects == "cor"){
    modMatrices$rho_randomEffects <- matrixsetup_rho(rho_randomEffects,
                                                     expcov=list(expRanEffects),
                                                     nNode = nCor,
                                                     nGroup = nGroup,
                                                     labels = corvars,
                                                     equal = FALSE,
                                                     sampletable = sampleStats,
                                                     name = "rho_randomEffects")
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

  # Extra matrices (both LVM and meta-analytic):
  model@extramatrices <- list(
    # LVM matrices (nNode dimension):
    D = psychonetrics::duplicationMatrix(nNode),
    L = psychonetrics::eliminationMatrix(nNode),
    Lstar = psychonetrics::eliminationMatrix(nNode, diag=FALSE),
    Dstar = psychonetrics::duplicationMatrix(nNode,diag = FALSE),
    In = as(diag(nNode),"dMatrix"),
    A = psychonetrics::diagonalizationMatrix(nNode),
    C = as(lavaan::lav_matrix_commutation(nNode, nLatent),"dMatrix"),
    C_chol = as(lavaan::lav_matrix_commutation(nNode, nNode),"dMatrix"),

    # LVM matrices (nLatent dimension):
    Deta = psychonetrics::duplicationMatrix(nLatent),
    L_eta = psychonetrics::eliminationMatrix(nLatent),
    Dstar_eta = psychonetrics::duplicationMatrix(nLatent,diag = FALSE),
    Inlatent = as(diag(nLatent),"dMatrix"),
    Cbeta = as(lavaan::lav_matrix_commutation(nLatent, nLatent),"dMatrix"),
    Aeta = psychonetrics::diagonalizationMatrix(nLatent),

    # Random effects matrices (nCor dimension):
    D_c = psychonetrics::duplicationMatrix(nCor),
    L_c = psychonetrics::eliminationMatrix(nCor),
    Lstar_c = psychonetrics::eliminationMatrix(nCor, diag=FALSE),
    Dstar_c = psychonetrics::duplicationMatrix(nCor,diag = FALSE),
    In_c = as(diag(nCor),"dMatrix"),
    A_c = psychonetrics::diagonalizationMatrix(nCor),
    C_c = as(lavaan::lav_matrix_commutation(nCor,nCor),"dMatrix"),

    # V matrices:
    V = avgVmat,
    Vall = Vmats,
    Vmethod = Vmethod,
    Vestimation = Vestimation
  )

  # Form the model matrices
  model@modelmatrices <- formModelMatrices(model)

  ### Baseline model ###
  if (baseline_saturated){
    # Dummy sample stats (not corinput for baseline/saturated):
    sampleStats2 <- sampleStats
    sampleStats2@corinput <- FALSE

    # Form baseline model:
    model@baseline_saturated$baseline <- varcov(data,
                                                mu = rep(0,nCor),
                                                type = "chol",
                                                lowertri = "diag",
                                                vars = corvars,
                                                missing = "pairwise",
                                                estimator = estimator,
                                                baseline_saturated = FALSE,
                                                sampleStats=sampleStats2)

    model@baseline_saturated$baseline@sample@fullFIML <- FALSE

    ### Saturated model ###
    model@baseline_saturated$saturated <- varcov(data,
                                                 type = "chol",
                                                 lowertri = "full",
                                                 vars = corvars,
                                                 missing = "pairwise",
                                                 estimator = estimator,
                                                 baseline_saturated = FALSE,
                                                 sampleStats=sampleStats2)

    model@baseline_saturated$saturated@sample@fullFIML <- FALSE

    # if not FIML, Treat as computed:
    if (estimator != "FIML"){
      model@baseline_saturated$saturated@computed <- TRUE
      model@baseline_saturated$saturated@objective <- psychonetrics_fitfunction(parVector(model@baseline_saturated$saturated),model@baseline_saturated$saturated)
    }
  }

  # Identify model (reuse LVM identification):
  model@identification <- identification
  model <- identify_lvm(model)

  if (missing(optimizer)){
    model <- setoptimizer(model, "default")
  } else {
    model <- setoptimizer(model, optimizer)
  }

  # FIXME: No C++ implementation yet, force R-only:
  model <- usecpp(model, FALSE)

  # Return model:
  return(model)
}
