
# Currently RI-CLPM is implemented using the lvm framework

ri_clpm <- function(
    data,
    vars,
    
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

  # Check stationary:
  # if (any(!stationary %in% c("none","contemporaneous","innovation","temporal","intercepts","random_intercept"))){
  #   stop("Stationary should include only 'none', 'random_intercept', 'contemporaneous', 'innovation', 'temporal', or 'intercepts'")
  # }
  
  # No NA in vars supported:
  if (any(is.na(vars))){
    stop("NA values in the design matrix ('vars') are not supported yet for the RI-CLPM.")
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
  
  # Submodel:
  mod@submodel <- "RI_CLPM"
  
  # Add saturated:
  mod@baseline_saturated$saturated <- varcov(data,
                                               type = "chol", 
                                               lowertri = "full", 
                                               vars = varsDF$column,
                                               ...,
                                               baseline_saturated = FALSE)
   
  # FIXME: Add baseline model
  
  # FIXME: Add constraints (flexible function to do this also in steps)
  
  # Return model:
  return(mod)
}
