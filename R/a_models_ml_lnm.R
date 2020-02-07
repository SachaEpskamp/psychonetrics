ml_lnm <- function(...){
  ml_lvm(...,within_latent = "ggm", between_latent = "ggm", within_residual = "cov", between_residual = "cov")
}