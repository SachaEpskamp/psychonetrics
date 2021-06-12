ml_ts_lvgvar <- function(...){
  ml_tsdlvm1(...,  within_latent = "ggm", between_latent = "ggm")
}

ml_var <- function(..., 
                   contemporaneous = c("cov","chol","prec","ggm"), 
                   between = c("cov","chol","prec","ggm")){
  
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  
  ml_tsdlvm1(...,  within_latent = contemporaneous, between_latent = between)
}

ml_gvar <- function(..., 
                   contemporaneous = c("ggm","cov","chol","prec"), 
                   between = c("ggm","cov","chol","prec")){
  
  contemporaneous <- match.arg(contemporaneous)
  between <- match.arg(between)
  
  ml_tsdlvm1(...,  within_latent = contemporaneous, between_latent = between)
}