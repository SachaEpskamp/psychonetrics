# ml_var / ml_gvar: the primary names of the multi-level VAR / graphical VAR
# family (ml_var1/ml_gvar1 are identical and remain fully supported).
#
# Historically (<= 0.15) these were wrappers around dlvm1() (full-information
# ML with every measurement occasion treated as a wave). As of 0.16.7 they are
# aliases for ml_var1() / ml_gvar1(), whose default estimator = "auto" keeps
# the 0.15 behavior for short (<= 10 occasions, panel-like) series - via the
# faster panelvar framework, with identical estimates - and switches to the
# two-level summary-statistics pseudo-ML for longer (ESM-like) series where
# full-information ML is infeasible.
#
# ml_ts_lvgvar (latent version) remains a deprecated wrapper for dlvm1().

ml_ts_lvgvar <- function(...){
  .Deprecated("dlvm1")
  dots <- list(...)
  if (!"estimator" %in% names(dots)) dots$estimator <- "FIML"
  do.call(dlvm1, c(dots, list(within_latent = "ggm", between_latent = "ggm")))
}

ml_var <- function(...,
                   contemporaneous = c("cov","chol","prec","ggm","cor"),
                   between = c("cov","chol","prec","ggm","cor")){
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  ml_var1(..., within_latent = contemporaneous, between_latent = between)
}

ml_gvar <- function(...,
                   contemporaneous = c("ggm","cov","chol","prec","cor"),
                   between = c("ggm","cov","chol","prec","cor")){
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  ml_var1(..., within_latent = contemporaneous, between_latent = between)
}
