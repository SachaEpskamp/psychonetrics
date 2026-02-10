# Deprecated wrappers: use dlvm1() or panellvgvar() directly instead.

ml_ts_lvgvar <- function(...){
  .Deprecated("panellvgvar")
  dots <- list(...)
  if (!"estimator" %in% names(dots)) dots$estimator <- "FIML"
  do.call(dlvm1, c(dots, list(within_latent = "ggm", between_latent = "ggm")))
}

ml_var <- function(...,
                   contemporaneous = c("cov","chol","prec","ggm"),
                   between = c("cov","chol","prec","ggm")){
  .Deprecated("dlvm1")
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  dots <- list(...)
  if (!"estimator" %in% names(dots)) dots$estimator <- "FIML"
  do.call(dlvm1, c(dots, list(within_latent = contemporaneous, between_latent = between)))
}

ml_gvar <- function(...,
                   contemporaneous = c("ggm","cov","chol","prec"),
                   between = c("ggm","cov","chol","prec")){
  .Deprecated("dlvm1")
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  dots <- list(...)
  if (!"estimator" %in% names(dots)) dots$estimator <- "FIML"
  do.call(dlvm1, c(dots, list(within_latent = contemporaneous, between_latent = between)))
}
