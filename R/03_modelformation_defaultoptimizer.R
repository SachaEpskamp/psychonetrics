# Single source of truth for the default optimizer. Model constructors,
# setoptimizer(..., "default") and runmodel's "default" branch all call this
# function, so that (for example) a refit()-ed model chooses the same
# optimizer as a freshly specified one. The argument x (a psychonetrics model)
# is currently unused but kept so model-dependent defaults can be added later.
defaultoptimizer <- function(x){
  # "cpp_L-BFGS-B"
  "nlminb"
}
