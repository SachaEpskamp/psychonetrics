update_baseline <- function(x, baseline){
  x@baseline_saturated$baseline <- baseline

  x
}

update_saturated <- function(x, saturated){
  x@baseline_saturated$saturated <- saturated

  x
}
