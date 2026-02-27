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

  # Determine number of variables and set up lambda/residuals:
  if (is.matrix(vars)){
    # Wide format: vars is a design matrix
    nVars <- nrow(vars)
    if (is.null(rownames(vars))){
      rownames(vars) <- paste0("Eta_",seq_len(nVars))
    }
    latentNames <- rownames(vars)
  } else if (is.character(vars)){
    # Long format: vars is a character vector of variable names
    nVars <- length(vars)
    latentNames <- vars
  } else {
    stop("'vars' must be a design matrix (wide format) or a character vector of variable names (long format).")
  }

  I <- diag(nVars)
  O <- matrix(0, nVars, nVars)

  # Conditionally pass data (allows covs-only usage without data):
  if (missing(data)) {
    dlvm1(vars = vars,
          lambda = I,
          within_latent = within_latent,
          within_residual = "cov", sigma_epsilon_within = O,
          between_latent = between_latent,
          between_residual = "cov", sigma_epsilon_between = O,
          latents = latentNames,
          ...
    )
  } else {
    dlvm1(data = data, vars = vars,
          lambda = I,
          within_latent = within_latent,
          within_residual = "cov", sigma_epsilon_within = O,
          between_latent = between_latent,
          between_residual = "cov", sigma_epsilon_between = O,
          latents = latentNames,
          ...
    )
  }
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

  # Determine number of variables and set up lambda/residuals:
  if (is.matrix(vars)){
    # Wide format: vars is a design matrix
    nVars <- nrow(vars)
    if (is.null(rownames(vars))){
      rownames(vars) <- paste0("Eta_",seq_len(nVars))
    }
    latentNames <- rownames(vars)
  } else if (is.character(vars)){
    # Long format: vars is a character vector of variable names
    nVars <- length(vars)
    latentNames <- vars
  } else {
    stop("'vars' must be a design matrix (wide format) or a character vector of variable names (long format).")
  }

  I <- diag(nVars)
  O <- matrix(0, nVars, nVars)

  # Conditionally pass data (allows covs-only usage without data):
  if (missing(data)) {
    dlvm1(vars = vars,
          lambda = I,
          within_latent = within_latent,
          within_residual = "cov", sigma_epsilon_within = O,
          between_latent = between_latent,
          between_residual = "cov", sigma_epsilon_between = O,
          latents = latentNames,
          ...
    )
  } else {
    dlvm1(data = data, vars = vars,
          lambda = I,
          within_latent = within_latent,
          within_residual = "cov", sigma_epsilon_within = O,
          between_latent = between_latent,
          between_residual = "cov", sigma_epsilon_between = O,
          latents = latentNames,
          ...
    )
  }
}


# Panel latent variable GVAR (new canonical name):
panellvgvar <- function(...){
  dlvm1(..., within_latent = "ggm", between_latent = "ggm")
}

# Deprecated alias (use panellvgvar instead):
panel_lvgvar <- function(...){
  .Deprecated("panellvgvar")
  panellvgvar(...)
}
