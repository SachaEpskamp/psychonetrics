matrixsetup_flexcov <- function(
  sigma,lowertri,omega,delta,kappa,
  type = "cov",
  name= "",
  sampleStats,
  equal = character(0),
  nNode,
  expCov,
  nGroup,
  labels
){
  modMatrices <- list()
  
  # type varcov:
  if (type == "cov"){
    mat <- paste0("sigma_",name)
    modMatrices[[mat]] <- matrixsetup_sigma(sigma, 
                                                name = mat,
                                                expcov=expCov,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = labels,
                                                equal = mat %in% equal, sampletable = sampleStats)    
  } else if (type == "chol"){
    mat <-  paste0("lowertri_",name)
    modMatrices[[mat]] <- matrixsetup_lowertri(lowertri, 
                                                      name = mat,
                                                      expcov=expCov,
                                                      nNode = nNode, 
                                                      nGroup = nGroup, 
                                                      labels = labels,
                                                      equal = mat %in% equal, sampletable = sampleStats)
  } else if (type == "ggm"){
    mat <- paste0("omega_",name)
    # Add omega matrix:
    modMatrices[[mat]] <- matrixsetup_omega(omega, 
                                                name = mat,
                                                expcov=expCov,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = labels,
                                                equal = mat  %in% equal, sampletable = sampleStats)
    
    # Add delta matrix:
    matDelta <- paste0("delta_",name)
    modMatrices[[matDelta]] <- matrixsetup_delta(delta, 
                                                name = matDelta,
                                                expcov=expCov,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = labels,
                                                equal = matDelta %in% equal, sampletable = sampleStats,
                                            omegaStart =  modMatrices[[mat]]$start) 
  } else if (type == "prec"){
    mat <- paste0("kappa_",name)
    # Add omega matrix:
    modMatrices[[mat]] <- matrixsetup_kappa(kappa, 
                                                name = mat,
                                                expcov=expCov,
                                                nNode = nNode, 
                                                nGroup = nGroup, 
                                                labels = labels,
                                                equal = mat %in% equal, sampletable = sampleStats)
  }
  modMatrices
}