panelvar <- function(data,vars,
                     within_latent = c("cov","chol","prec","ggm"), 
                     between_latent = c("cov","chol","prec","ggm"), 
                     ...){
  # Match arg:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  
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
        within_latent = within_latent,
        within_residual = "cov", sigma_epsilon_within = O,
        between_latent = between_latent,
        between_residual = "cov", sigma_epsilon_between = O,
        latents = rownames(vars),
        ...
        )
}

panelgvar <- function(data,vars,
                      within_latent = c("ggm","chol","cov","prec"), 
                      between_latent = c("ggm","chol","cov","prec"), 
                      ...){
  # Match arg:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)
  
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
        within_latent = within_latent,
        within_residual = "cov", sigma_epsilon_within = O,
        between_latent = between_latent,
        between_residual = "cov", sigma_epsilon_between = O,
        latents = rownames(vars),
        ...
  )
}


# Panel gvar:
panel_lvgvar <- function(...){
  dlvm1(..., within_latent = "ggm", between_latent = "ggm")
}

