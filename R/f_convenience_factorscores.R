factorscores <- function(data, model,
      method = c('bartlett','regression')){

  method <- match.arg(method)
  
  # Select variables:
  if (!all(model@sample@variables$label %in% names(data))){
    stop("Not all variables of model are present in data")
  }
  data <- data[,model@sample@variables$label]
    
  # Check input:
  if (any(is.na(data))){
    stop("Missing data is not supported at the moment.")
  }
  
  if (!is(model,"psychonetrics")){
    stop("Input must be a 'psychonetrics' object.")
  }
  if (model@model != "lvm"){
    stop("Only 'lvm' framework is currently supported.")
  }
  
  # Number of groups:
  nGroups <- nrow(model@sample@groups)
  
  if (nGroups > 1){
    stop("Currently only single group models are supported.")
  }
  
  # Extract matrices:
  if (model@model == "lvm"){
    
    g <- 1
    
    sigma_zeta <- model@modelmatrices[[g]]$sigma_zeta
    beta <- model@modelmatrices[[g]]$beta
    sigma_epsilon <- model@modelmatrices[[g]]$sigma_epsilon
    lambda <- model@modelmatrices[[g]]$lambda
    sigma <- model@modelmatrices[[g]]$sigma
    kappa <- model@modelmatrices[[g]]$kappa
    mu <- model@modelmatrices[[g]]$mu
    
    # Weights:
    if (method == "regression"){
      
      W <- sigma_zeta %*% t(lambda) %*% kappa
      
    } else if (method == "bartlett") {
      
      W <- solve_symmetric_cpp_matrixonly(t(lambda) %*% solve_symmetric_cpp_matrixonly(sigma_epsilon) %*% lambda) %*% t(lambda) %*% solve_symmetric_cpp_matrixonly(sigma_epsilon)
      
    }
    
  }
  
  # PRedict:
  eta <- t(W %*% (t(data) - as.vector(mu)))
  
  return(eta)
  
}