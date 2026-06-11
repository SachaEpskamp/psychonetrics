# Simple wrapper around lvm
bifactor <- function(data, lambda, latents,  bifactor = "g", ...){
  if (missing (lambda)){
    stop("lambda argument may not be missing")
  }
  
  # Number of latents:
  nLatent <- ncol(lambda)
  
  
  # If latents is not provided, make it:
  if (missing(latents)){
    latents <- paste0("Eta_",seq_len(nLatent))
  }
  if (length(latents) != nLatent){
    stop("Length of 'latents' is not equal to number of latent variables in model.")
  }
  
  # Augment lambda:
  lambda <- cbind(lambda,1)
  
  latents <- c(latents, bifactor)
  
  # Return model:
  mod <- lvm(data, lambda = lambda,  latents=latents, ...,
      sigma_zeta = "diag", #
      kappa_zeta = "diag", # Precision
      omega_zeta = "zero", # Partial correlations
      lowertri_zeta = "diag", # Cholesky
      delta_zeta = "diag")
  
  # Deterministic starting values (previously non-reproducible rnorm draws).
  # Small non-zero starts break the symmetry between the general and group
  # factors without depending on the RNG state.
  # Adjust factor loading start values:
  ind <- mod@parameters$matrix == "lambda" & !mod@parameters$fixed
  mod@parameters$est[ind] <- 0.01

  # Adjust other start values:
  ind <- mod@parameters$matrix != "lambda" &mod@parameters$matrix != "nu" & !mod@parameters$fixed & mod@parameters$row != mod@parameters$col
  mod@parameters$est[ind] <-  0.01
  ind <- mod@parameters$matrix != "lambda" &mod@parameters$matrix != "nu" & !mod@parameters$fixed & mod@parameters$row == mod@parameters$col
  mod@parameters$est[ind] <- 0.5
  
  
  # Return:
  mod
}