# Fit function per group:
fimlEstimator_Gauss_group <- function(mu,sigma,data,kappa,...){
  curF <- 0 
  

  
  n <- nrow(data)
  if (any(eigen(kappa)$values < 0)){
    kappa <- Matrix::nearPD(kappa)$mat
  }
  logdetK <- max(log(det(kappa)),0)
  
  # For every subject:
   for (i in seq_len(n)){
     if (!any(is.na(data[i,]))){
       kappa_p <- kappa
       mu_p <- mu
       y <- unlist(data[i,])
       logdetK_p <- logdetK
     } else {
       obs <- unlist(!is.na(data[i,]))
       if (!any(obs)) next
       # Handle possible non positive definiteness:
       sig_p <- as.matrix(sigma)[obs,obs]
       if (any(eigen(sig_p)$values < 0)){
         logdetK_p <- 0
         if (all(eigen(sig_p)$values < 0)){
           next # FIXME: Not sure what to do else about this...
         }
         sig_p <- Matrix::nearPD(sig_p)$mat
         kappa_p <- solve(sig_p)
       } else {
         kappa_p <- solve(sig_p)
         logdetK_p <-max(log(det(kappa_p)),0)
       }

       
       mu_p <- mu[obs,]
       y <- unlist(data[i,obs])
     }
     if (any(eigen(kappa_p)$values < 0)){
       kappa_p <- Matrix::nearPD(kappa_p)$mat
     }
     curF <- curF + 1/n * (sum(diag(kappa_p %*% (y - mu_p) %*% t(y - mu_p)))  - logdetK_p)
   }
  
  as.numeric(curF)
}

# Fit function for Gauss ML: -2n* log likelihood
fimlEstimator_Gauss <- function(x, model){
  # Prepare
  prep <- prepareModel(x, model)


  # Fit function per group:
  fit_per_group <- prep$nPerGroup / prep$nTotal * sapply(prep$groupModels,do.call,what=fimlEstimator_Gauss_group)

  # Sum and return:
  sum(fit_per_group)
}