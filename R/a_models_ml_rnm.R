ml_rnm <- function(...){
  ml_lvm(..., within_latent = "cov", between_latent = "cov", within_residual = "ggm", between_residual = "ggm")
}