emergencystart <- function(x){
  
  ## For the ggm model, we can use glasso to help here! Only for single group now...
  nGroups <- nrow(x@sample@groups)
  
  if (x@model == 'varcov' && x@types[['y']] == 'ggm' && !is.null(x@baseline_saturated$saturated) && nGroups == 1){
    # Form matrices:
    mats <- formModelMatrices(x)
    
    # For each group, compute glasso and set starting values:
    nGroups <- nrow(x@sample@groups)
    # for (g in seq_len(nGroups)){
    g <- 1 # FIXME: Add for multiple groups!
      
      # Extract saturated covs:
      satCovs <- spectralshift(x@baseline_saturated$saturated@modelmatrices[[g]]$sigma)
      
      # Model zeroes:
      net <- mats[[g]]$omega != 0
      zeroes <- which(net & diag(ncol(net)) != 1, arr.ind = TRUE)
      
      # Glasso result:
      if (nrow(zeroes) > 0){
        suppressWarnings(glas <- glasso::glasso(as.matrix(satCovs), 0, zero = zeroes))     
        omega <- -1*cov2cor(glas$wi)
        diag(omega) <- 0
      } else {
        omega <- qgraph::wi2net(solve(as.matrix(satCovs)))
      }

      
      # To omega:
    
      
      # Overwrite start:
      for (i in which(x@parameters$matrix == "omega" & x@parameters$group_id == g & !x@parameters$fixed)){
        x@parameters$est[i] <- omega[x@parameters$row[i], x@parameters$col[i]]
      }
    # }
      # Start for delta:
      if (!x@sample@corinput){
        # Scale to correlations:
        toR <- diag(1/sqrt(diag(solve_symmetric(diag(nrow(omega)) - omega))))
        
        # Scale to covariances:
        toCov <- diag(sqrt(diag(satCovs)))
        
        # This makes delta:
        delta <- toCov %*% toR
        
        # Now add these:
        # Overwrite start:
        for (i in which(x@parameters$matrix == "delta" & x@parameters$group_id == g & !x@parameters$fixed)){
          x@parameters$est[i] <- delta[x@parameters$row[i], x@parameters$col[i]]
        }
      }
  } else {
    # Try fixing:
    # For all matrices in the model:
    allMats <- x@matrices[!grepl("(mu)|(nu)",x@matrices$name),]
    for (mat in seq_len(nrow(allMats))){
      # Adjust starting values. Identity matrix for posdef matrices, Near-null matrix for anything else:
      if (allMats$posdef[mat]){
        x@parameters$est[x@parameters$matrix == allMats$name[mat] & !x@parameters$fixed] <- 
          ifelse(x@parameters$row[x@parameters$matrix == allMats$name[mat] & !x@parameters$fixed] ==
                   x@parameters$col[x@parameters$matrix == allMats$name[mat] & !x@parameters$fixed], 1, 0)
      } else {
        x@parameters$est[x@parameters$matrix == allMats$name[mat]& !x@parameters$fixed] <- 
          (x@parameters$est[x@parameters$matrix == allMats$name[mat]& !x@parameters$fixed] != 0) * 1e-7 * sign(x@parameters$est[x@parameters$matrix == allMats$name[mat]& !x@parameters$fixed])
      }
      
    }
    
  }
  
 
  
  return(x)
}