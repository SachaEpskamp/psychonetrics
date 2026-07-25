matrixsetup_mu <- function(
  mu, # sigma argument
  nNode, # Number of nodes
  nGroup, # Number of groups
  expmeans, # Expected covariance matrices (list).
  labels,
  equal = FALSE,
  sampletable,
  name = "mu",
  meanstructure = TRUE,
  # Wave-specific ("free") mean structures (means = "free" in the panelvar and
  # dlvm1 frameworks). With nTime = 1 (the default) the matrix is the usual
  # nNode * nGroup one and everything below is exactly as it always was; with
  # nTime > 1 the matrix becomes nNode * nTime * nGroup, with one (free) mean
  # per variable--wave combination present in the design:
  nTime = 1,
  design, # nNode * nTime design matrix (1 = variable observed at that wave)
  timelabels # Column labels used for the wave-specific case
){

  # Stationary (single mean per variable) case:
  if (nTime == 1){
    # Fix mu
    mu <- fixMu(mu,nGroup,nNode,equal)

    # For each group, form starting values:
    muStart <- mu
    for (g in 1:nGroup){
      # Current estimate:
      meanest <-  as.matrix(expmeans[[g]])

      # Covs with sample covs:
      muStart[,g] <- 1*(mu[,g]!=0) * meanest
    }

    # Form the model matrix part:
    return(list(mu,
         mat =  name,
         op =  "~1",
         symmetrical= FALSE,
         sampletable=sampletable,
         rownames = labels,
         colnames = "1",
         start = muStart,
         allFixed = !meanstructure))
  }

  ### Wave-specific means ###
  if (missing(design)) design <- matrix(1, nNode, nTime)
  if (missing(timelabels)) timelabels <- paste0("T",seq_len(nTime))

  # Default pattern: one free mean per variable--wave combination in the design:
  if (missing(mu) || is.null(mu)){
    mu <- 1*(design == 1)
  } else if (is.character(mu) && length(mu) == 1){
    mu <- switch(mu,
                 "full" = 1*(design == 1),
                 "free" = 1*(design == 1),
                 "zero" = matrix(0,nNode,nTime),
                 "empty" = matrix(0,nNode,nTime),
                 stop(paste0("'",name,"' can only be 'full', 'zero' or a model matrix"))
    )
  } else if (is.vector(mu) && length(mu) == nNode){
    # A single pattern per variable: repeat over waves:
    mu <- matrix(mu, nNode, nTime)
  }

  # Fix the matrix (handles matrix / array input and equality across groups):
  mu <- fixMatrix(mu, nGroup = nGroup, nrows = nNode, ncols = nTime, equal = equal)

  # Variable--wave combinations that are not in the design have no mean to
  # estimate:
  if (any(design == 0)){
    for (g in seq_len(nGroup)){
      cur <- mu[,,g]
      if (any(cur[design == 0] != 0)){
        warning(paste0("Elements of '",name,"' at variable--wave combinations that are not in the design matrix are fixed to zero."))
        cur[design == 0] <- 0
        mu[,,g] <- cur
      }
    }
  }

  # For each group, form starting values:
  muStart <- mu
  for (g in 1:nGroup){
    # Current estimate (a nNode * nTime matrix, or a vector to be repeated):
    meanest <- as.matrix(expmeans[[g]])
    if (ncol(meanest) == 1 && nTime > 1){
      meanest <- matrix(meanest[,1], nNode, nTime)
    }
    meanest[!is.finite(meanest)] <- 0

    muStart[,,g] <- 1*(mu[,,g]!=0) * meanest
  }

  # Form the model matrix part:
  list(mu,
       mat =  name,
       op =  "~1",
       symmetrical= FALSE,
       sampletable=sampletable,
       rownames = labels,
       colnames = timelabels,
       start = muStart,
       allFixed = !meanstructure)
}
