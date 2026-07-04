# Panel graphical VAR: the panelvar model with GGM parameterizations for the
# within- and between-person structures. Since psychonetrics 0.16.2 this uses
# the dedicated panelvar framework (see a_models_panelvar.R) rather than
# dlvm1 with dummy matrices:
panelgvar <- function(data,vars,
                      within_latent = c("ggm","chol","cov","prec","cor"),
                      between_latent = c("ggm","chol","cov","prec","cor"),
                      ...){
  # Match arg:
  within_latent <- match.arg(within_latent)
  between_latent <- match.arg(between_latent)

  # Conditionally pass data (allows covs-only usage without data):
  if (missing(data)) {
    panelvar(vars = vars,
             within_latent = within_latent,
             between_latent = between_latent,
             ...)
  } else {
    panelvar(data = data, vars = vars,
             within_latent = within_latent,
             between_latent = between_latent,
             ...)
  }
}


# Panel latent variable GVAR (latent variables: uses the dlvm1 framework):
panellvgvar <- function(...){
  dlvm1(..., within_latent = "ggm", between_latent = "ggm")
}

# Deprecated alias (use panellvgvar instead):
panel_lvgvar <- function(...){
  .Deprecated("panellvgvar")
  panellvgvar(...)
}
