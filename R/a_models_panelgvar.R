panelvar <- function(data,vars,...){
  # Extract var names:
  if (missing(vars)){
    stop("'vars' argument may not be missing")
  }
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix, with rows indicating variables and columns indicating measurements.")
  }
  
  
  I <- diag(nrow(vars))
  O <- matrix(0,nrow(vars),nrow(vars))
  
  
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("Eta_",seq_len(nrow(vars)))
  }
  
  dlvm1(data,vars,
        lambda = I,
        within_latent = "chol",
        within_residual = "cov", sigma_epsilon_within = O,
        between_latent = "chol",
        between_residual = "cov", sigma_epsilon_between = O,
        latents = rownames(vars),
        ...
        )
}

panelgvar <- function(data,vars,...){
  # Extract var names:
  if (missing(vars)){
    stop("'vars' argument may not be missing")
  }
  if (!is.matrix(vars)){
    stop("'vars' must be a design matrix, with rows indicating variables and columns indicating measurements.")
  }
  
  
  I <- diag(nrow(vars))
  O <- matrix(0,nrow(vars),nrow(vars))
  
  if (is.null(rownames(vars))){
    rownames(vars) <- paste0("Eta_",seq_len(nrow(vars)))
  }
  
  dlvm1(data,vars,
        lambda = I,
        within_latent = "ggm",
        within_residual = "cov", sigma_epsilon_within = O,
        between_latent = "ggm",
        between_residual = "cov", sigma_epsilon_between = O,
        latents = rownames(vars),
        ...
  )
}


# Panel gvar:
panel_lvgvar <- function(...){
  dlvm1(..., within_latent = "ggm", between_latent = "ggm")
}

