# Deprecated: use dlvm1() with long-format data instead.
# This function is retained for backward compatibility.
ml_tsdlvm1 <- function(
  data,
  beepvar,
  idvar,
  vars,
  groups,
  estimator = "FIML",
  standardize = c("none","z","quantile"),
  ...
){
  .Deprecated("dlvm1")
  standardize <- match.arg(standardize)

  dlvm1(data = data, vars = vars, idvar = idvar, beepvar = beepvar,
        groups = groups, estimator = estimator, standardize = standardize,
        datatype = "long", ...)
}
