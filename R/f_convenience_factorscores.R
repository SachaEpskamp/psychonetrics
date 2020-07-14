factorscores <- function(data, model,
                         method = c('bartlett','regression')){
  
  method <- match.arg(method)
  
  # Number of groups:
  nGroups <- nrow(model@sample@groups)
  
  if (nGroups == 1){
    data[['GROUPDUMMY']] <- "fullsample"
    groupvar <- "GROUPDUMMY"
  } else {
    groupvar <- model@sample@groupvar
  }
  
  
  # Select variables:
  if (!all(model@sample@variables$label %in% names(data))){
    stop("Not all variables of model are present in data")
  }
  obsVars <- model@sample@variables$label
  data <- data[,c(obsVars,groupvar)]
  
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
  
 
  
  # Latent names:
  latNames <- unique(model@parameters$var2[model@parameters$matrix=="lambda"])
  
  # Dummy matrix with responses:
  eta <- matrix(NA,nrow(data),length(latNames))
  
  # Extract matrices:
  if (model@model == "lvm"){
    for (g in seq_len(nrow(model@sample@groups))){
      
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
      
      # PRedict:
      groupID <- model@sample@groups$label[g]
      eta[data[[groupvar]] == groupID,] <- t(W %*% (t( data[data[[groupvar]] == groupID,obsVars]) - as.vector(mu)))
    }
  } else {
    stop("Only 'lvm' models are currently supported.")
  }
  
  eta <- as.data.frame(eta)
  names(eta) <- latNames
  if (nGroups > 1){
    eta[[groupvar]] <- data[[groupvar]]
  }
  
  return(eta)
  
}