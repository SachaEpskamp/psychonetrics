# Implied model for precision. Requires appropriate model matrices:
implied_Ising <- function(model, all = FALSE){
  if (model@cpp){
    x <- formModelMatrices_cpp(model)
  } else {
    x <- formModelMatrices(model)  
  }
  
  
    # For each group:
    nGroup <- length(x)
    
    for (g in 1:nGroup){
      if (model@types$beta_model == "log_beta"){
        x[[g]]$beta <- exp(x[[g]]$log_beta)
      } else {
        x[[g]]$log_beta <- log(x[[g]]$beta)
      }
    }

    # When all = TRUE, compute model-implied means and covariances
    # from the Ising expectations (used after optimization for fit measures):
    if (all) {
      for (g in 1:nGroup) {
        exp <- isingExpectation(
          x[[g]]$omega, x[[g]]$tau, x[[g]]$beta,
          model@extramatrices$responses,
          model@extramatrices$min_sum
        )
        x[[g]]$mu <- exp$exp_v1
        x[[g]]$sigma <- exp$exp_v2 - tcrossprod(exp$exp_v1)
      }
    }

  x
}
