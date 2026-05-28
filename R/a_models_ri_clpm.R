
# Currently RI-CLPM is implemented using the lvm framework

ri_clpm <- function(
    data,
    vars,

    # Data format (mirrors dlvm1):
    datatype = c("auto", "wide", "long"), # auto-detected from 'vars'
    idvar,    # Subject ID variable (required for long format)
    beepvar,  # Measurement occasion variable (optional for long format)

    # Standardization:
    standardize = c("none", "z", "quantile", "z_per_wave"),

    # Lambda (not yet supported):
    lambda,

    # Types:
    type = c("cov","chol","prec","ggm"),

    # Which should be stationary:
    # stationary = "random_intercept", # "none", random_intercept", "contemporaneous", "innovation", "temporal" "intercepts" accepted

    # Verbose:
    verbose = FALSE, # Verbose output

    ... # Arguments used in lvm(...)
    ){

  # Types to use:
  type <- match.arg(type)
  datatype <- match.arg(datatype)
  standardize <- match.arg(standardize)

  # Track whether data was supplied:
  .data_missing <- missing(data)

  # CRAN-check workarounds for long-format conversion:
  variable <- NULL
  value <- NULL

  # --- Auto-detect data format (like dlvm1) ---
  if (missing(vars)){
    # If vars is missing but idvar is provided, assume long format:
    if (!missing(idvar)){
      datatype <- "long"
      if (is.matrix(data)) data <- as.data.frame(data)
      if (!is.data.frame(data)) stop("'data' must be a data frame")
      vars <- names(data)
      vars <- vars[vars != idvar]
      if (!missing(beepvar)) vars <- vars[vars != beepvar]
    } else {
      stop("'vars' argument may not be missing")
    }
  }
  if (datatype == "auto"){
    if (is.matrix(vars)){
      datatype <- "wide"
    } else if (is.character(vars)){
      datatype <- "long"
    } else {
      stop("Could not auto-detect data type from 'vars'. Please specify datatype = 'wide' or 'long'.")
    }
  }

  # --- Long-to-wide conversion (mirrors dlvm1) ---
  if (datatype == "long"){
    if (missing(idvar)){
      stop("'idvar' is required when datatype = 'long' or when 'vars' is a character vector.")
    }
    if (is.matrix(data)) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("'data' must be a data frame")
    if (!is.character(vars) || is.matrix(vars)){
      stop("When datatype = 'long', 'vars' must be a character vector of variable names.")
    }
    if (!idvar %in% names(data)){
      stop("'idvar' argument does not correspond to column name of 'data'")
    }

    # Remove rows with NA cluster:
    if (any(is.na(data[[idvar]]))){
      warning("Rows with NA in idvar removed.")
      data <- data[!is.na(data[[idvar]]),]
    }

    # Handle beepvar:
    if (missing(beepvar)){
      data <- data[order(data[[idvar]]),]
      beepvar <- "BEEPVAR"
      data[[beepvar]] <- unlist(tapply(data[[idvar]], data[[idvar]], seq_along))
    }
    if (any(data[[beepvar]][!is.na(data[[beepvar]])] %% 1 != 0)){
      stop("'beepvar' does not encode integer values.")
    }
    if (any(tapply(data[[beepvar]], data[[idvar]], function(x) any(duplicated(x)), simplify = TRUE))){
      stop("'beepvar' contains duplicate values for one or more cases.")
    }

    maxInCluster <- max(data[[beepvar]])

    # Reshape long to wide:
    datalong <- tidyr::gather(data, variable, value, vars)
    datawide <- tidyr::pivot_wider(datalong, id_cols = idvar,
                                   values_from = "value", names_from = c("variable", beepvar))

    # Construct design matrix:
    rowVars <- vars
    colVars <- as.character(seq(maxInCluster))
    designLong <- matrix(NA, length(rowVars), length(colVars))
    for (i in seq_along(rowVars)){
      for (j in seq_along(colVars)){
        varName <- paste0("^", rowVars[i], "_", colVars[j], "$")
        if (length(which(grepl(varName, names(datawide)))) == 1){
          designLong[i, j] <- paste0(rowVars[i], "_", colVars[j])
        }
      }
    }
    rownames(designLong) <- rowVars

    data <- datawide
    vars <- designLong
  }

  # --- Validate vars is now a design matrix ---
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix (rows = variables, columns = measurements) or a character vector of variable names when using long-format data.")
  }

  # --- Standardization (per variable across waves, or per variable x wave) ---
  if (standardize != "none" && !.data_missing && !is.null(data)){
    # Coerce to a plain data.frame so single-column extraction returns vectors
    # (tibbles from pivot_wider would otherwise return 1-column tibbles):
    data <- as.data.frame(data)
    for (v in seq_len(nrow(vars))){
      varCols <- na.omit(vars[v, ])
      if (length(varCols) == 0) next
      if (standardize == "z"){
        # Pool all waves of this variable, single mean/SD:
        allValues <- unlist(data[, varCols, drop = FALSE])
        m <- mean(allValues, na.rm = TRUE); s <- sd(allValues, na.rm = TRUE)
        for (col in varCols) data[[col]] <- (data[[col]] - m) / s
      } else if (standardize == "z_per_wave"){
        # Standardize each variable x wave column separately:
        for (col in varCols) data[[col]] <- (data[[col]] - mean(data[[col]], na.rm = TRUE)) / sd(data[[col]], na.rm = TRUE)
      } else if (standardize == "quantile"){
        allValues <- unlist(data[, varCols, drop = FALSE])
        allTransformed <- quantiletransform(allValues)
        n <- nrow(data)
        for (ci in seq_along(varCols)) data[[varCols[ci]]] <- allTransformed[((ci - 1) * n + 1):(ci * n)]
      }
    }
  }

  # Structural missing waves not yet supported (NA in the design matrix):
  if (any(is.na(vars))){
    stop("Structural missing waves (NA in the design matrix 'vars') are not yet supported for the RI-CLPM. Each variable must be measured at every wave. (Incomplete cases / dropout are handled automatically via FIML.)")
  }

  # Design matrix:
  # design <- as(1*(!is.na(vars)),"dMatrix")
  design <- as.matrix(1*(!is.na(vars)))
  vars <- as.matrix(vars)
  
  # Check time points:
  if (ncol(vars) < 3){
    stop("At least three time points are needed for the CI-RLPM")
  }
  
  # time per var:
  timePerVar <- as.vector(design * row(design))
  timePerVar <- timePerVar[timePerVar!=0]
  
  # Number of variables:
  nVar <- nrow(vars)
  
  # Number of measurements:
  nTime <- ncol(vars)
  
  # row names:
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("V",seq_len(nrow(vars)))
  }
  varnames <- rownames(vars)
  
  # col names:
  if (is.null(colnames(vars))){
    colnames(vars) <- paste0("T",seq_len(ncol(vars)))
  }
  timenames <- colnames(vars)
  
  # Data frame of the variables::
  varsDF <- data.frame(
    column = c(as.matrix(vars)),
    variable = rep(varnames,nTime),
    time = rep(timenames,each=nVar)
  )
    
    
  # Check lambda:
  if (missing(lambda)){
    if (verbose){
      message("'lambda' is missing, creating observed data only model.")
    }
    
    # FIXME: Not sure if needed
    
    # lambda <- diag(nVar)
    # O <- matrix(0, nVar, nVar)
    # omega_epsilon_within <- O
    # delta_epsilon_within <- O
    # kappa_epsilon_within <- O
    # sigma_epsilon_within <- O
    # lowertri_epsilon_within <- O
    # 
    # omega_epsilon_between <- O
    # delta_epsilon_between <- O
    # kappa_epsilon_between <- O
    # sigma_epsilon_between <- O
    # lowertri_epsilon_between <- O
  } else {
    stop("Latent variables are not yet supported for the RI-CLPM")
  }
  
  # The latents consist of the within-person deviations of each variable, here labeled C for contemporaneous:
  innovations <- paste0("C_",varsDF$column)
  n_observed <- nrow(varsDF)
  
  # And Random intercepts:
  RIs <- paste0("RI_",varnames)
  
  # all latents are:
  latents <- c(innovations,RIs)
  n_latent <- length(latents)
  
  # The beta matrix should be of this dimension
  beta <- matrix(0,n_latent,n_latent)
  
  # FIXME: Add equality constraints in model or after? Probably better after for better starting values?
  
  # # Random intercepts:
  # for (i in seq_len(nTime)){
  #   beta[(i-1)*nVar + seq_len(nVar),n_observed + seq_len(nVar)] <- diag(nVar)
  # }
  
  # Temporal effects:
  for (i in seq_len(nTime-1)){
    beta[i*nVar + seq_len(nVar),(i-1)*nVar + seq_len(nVar)] <- 1
  }
 
  # Contemporaneous effects:
  sigma_zeta <- diag(1,n_latent)
  
  # Fill covariances: + 1 for RI
  for (i in seq_len(nTime+1)){
    sigma_zeta[(i-1)*nVar + seq_len(nVar),(i-1)*nVar + seq_len(nVar)] <- 1
  }
  
  # Factor loadings:
  lambda <- matrix(0,n_observed,n_latent)
  
  # Fill RIs:
  for (i in seq_len(nTime)){
    lambda[(i-1)*nVar + seq_len(nVar),n_observed+seq_len(nVar)] <- diag(nVar)
  }
  
  # Fill contemporaneous:
  lambda[seq_len(n_observed),seq_len(n_observed)] <- diag(n_observed)
  
  
  # Form the base model:
  mod <- lvm(
    data = data,
    vars = varsDF$column,
    latent = type,
    
    # Latent names:
    latents = latents,
    identify = FALSE,
    
    # Lambda:
    lambda,
    
    
    # Structural:
    beta = beta,
    
    # Variances:
    sigma_zeta = sigma_zeta,
    omega_zeta = sigma_zeta,
    kappa_zeta = sigma_zeta,
    lowertri_zeta = ifelse(col(sigma_zeta) > row(sigma_zeta), 0, sigma_zeta),
    
    # Residual:
    residual = "cov",
    sigma_epsilon = "zero",
    
    # No baseline and saturated:
    baseline_saturated = FALSE,
    
    ...
    
  )
  
  # Fix all factorloadings to 1:
  which_fix <- mod@parameters$matrix=="lambda" & !mod@parameters$fixed
  mod@parameters$est[which_fix] <- 1
  mod@parameters$par[mod@parameters$matrix=="lambda"] <- 0
  mod@parameters$fixed[mod@parameters$matrix=="lambda"] <- TRUE
  mod@parameters$identified[mod@parameters$matrix=="lambda"] <- TRUE
  mod@parameters <- clearpars(mod@parameters,which_fix)
  
  # Identify nu_eta:
  mod@parameters$est[mod@parameters$matrix=="nu_eta"] <- 0
  mod@parameters$par[mod@parameters$matrix=="nu_eta"] <- 0
  mod@parameters$fixed[mod@parameters$matrix=="nu_eta"] <- TRUE
  mod@parameters$identified[mod@parameters$matrix=="nu_eta"] <- TRUE
  mod@parameters <- clearpars(mod@parameters,mod@parameters$matrix=="nu_eta")
  
  # Relabel:
  mod@parameters <- parRelabel(mod@parameters)

  
  # Simple start values for now:
  lt <- lower.tri(matrix(0,n_observed,n_observed))
  
  
  # Add extra matrices:
  mod@extramatrices$vars <- vars
  mod@extramatrices$varsDF <- varsDF
  mod@extramatrices$ri_clpm_type <- type
  
  # Submodel:
  mod@submodel <- "RI_CLPM"
  
  # Add saturated (all observed variances/covariances free):
  mod@baseline_saturated$saturated <- varcov(data,
                                               type = "chol",
                                               lowertri = "full",
                                               vars = varsDF$column,
                                               ...,
                                               baseline_saturated = FALSE)

  # Add baseline (independence: free observed variances and means, no covariances):
  mod@baseline_saturated$baseline <- varcov(data,
                                            type = "chol",
                                            lowertri = "diag",
                                            vars = varsDF$column,
                                            ...,
                                            baseline_saturated = FALSE)
  
  # Return model:
  return(mod)
}
