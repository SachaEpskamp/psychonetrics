latentgrowth <- function(
  vars,  # design matrix
  time = seq_len(ncol(vars)) - 1,
  covariates = character(0),
  covariates_as = c("regression","covariance"),
  ...){
  covariates_as <- match.arg(covariates_as)
  if (missing(vars)){
    stop("'vars' argument may not be missing")
  }
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix, with rows indicating variables and columns indicating measurements.")
  }
  
  
  # If rownames are missing, add them:
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("V",seq_len(nrow(vars)))
  }
  varnames <- rownames(vars)
  allVars <- c(na.omit(c(t(vars))),covariates)
  
  # Number of covariates:
  nCovariate <- length(covariates)
  
  # Number of regular variables:
  nReg <- length(allVars) - nCovariate
  
  # Number of intercepts and slopes:
  nIntSlope <- nrow(vars) * 2
  
  # Latent names:
  latNames <- c(paste0("int_",varnames),paste0("slope_",varnames),covariates)
  
  # Construct lambda:
  Lambda <- matrix(0, length(allVars), length(latNames))
  
  # Fill int and slopes:
  # Lambda[seq_len(nReg),seq_len(nIntSlope)] <- 1
  
  # Covariates?
  if (length(covariates) > 0){
    for (i in seq_along(covariates)){
      Lambda[allVars == covariates[i],latNames == covariates[i]] <- 1
    }
  }
  
  # Constrain factor loadings:
  for (i in 1:nrow(vars)){
    Lambda[allVars %in% na.omit(vars[i,]),i] <- 1
    Lambda[allVars %in% na.omit(vars[i,]),i+nrow(vars)] <- 1
  }

  # If the covariates are regressions, add these:
  nLat <- length(latNames)
  nLat_growth <- 2*length(varnames)
  beta <- matrix(0,nLat,nLat)
  if (covariates_as == "regression"){
    beta[seq_len(nLat_growth),-seq_len(nLat_growth)] <- 1

    psi <- diag(nLat)
    psi[seq_len(nLat_growth),seq_len(nLat_growth)] <- 1
  } else {
    psi <- matrix(1,nLat, nLat)
  }
    
  # Form model:
  mod <- lvm(lambda=Lambda, vars = allVars, latents = latNames, 
             identify = FALSE, nu = rep(0,length(allVars)),
             simplelambdastart = TRUE,beta=beta,
             sigma_zeta = psi,omega_zeta = psi,
             lowertri_zeta = psi, kappa_zeta = psi,
             ...)

  
  # Constrain factor loadings:
  for (i in 1:nrow(vars)){
    mod <- mod %>% fixpar("lambda",row=na.omit(vars[i,]),col=i,value = 1, verbose = FALSE)  
    for (j in seq_along(time)){
      if (!is.na(time[j])){
        if (!is.na(vars[i,j])){
          mod <- mod %>% fixpar("lambda",row=na.omit(vars[i,j]),col=i+nrow(vars),value = time[j], verbose = FALSE)   
        }
      }  
    }
  }
  


  
  # Identify:
  mod <- identify_lvm(mod)

  
  # FIXME: Silly start values:
  # mod <- emergencystart(mod)
  
  mod
}
