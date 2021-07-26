# Multi-level (random intercept) tsdlvm1 specification. This is the same as dlvm1(), but with the mlVAR like specification:
ml_tsdlvm1 <- function(
  data, 
  beepvar,
  # dayvar,
  idvar,
  vars, 
  groups,
  estimator = "FIML",
  standardize = c("none","z","quantile"),
  # nightskip = 1,
  ...
){
  # CRAN Check workarounds (sorry):
  . <- NULL
  variable <- NULL
  value <- NULL
  
  standardize <- match.arg(standardize)
  
  if (estimator != "FIML") stop("Only 'FIML' supported currently.")
  
  # Check idvar:
  if (missing(idvar)){
    stop("'idvar' may not be missing (use tsdlvm1 family instead).")
  }
  
  # Check if data is data frame:
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  
  # Check vars:
  if (missing(vars)){
    if (is.null(names(data))){
      stop("Dataset contains no column names.")
    }
    vars <- names(data)
    vars <- vars[vars!=idvar]
  } else {
    if (is.null(names(data))){
      stop("Dataset contains no column names.")
    }
  }
  
  # Number of idvar:
  if (!idvar %in% names(data)){
    stop("'idvar' argument does not correspond to column name of 'data'")
  }
  
  
  
  # Remove data with NA cluster:
  if (any(is.na(data[[idvar]]))){
    warning("Rows with NA cluster removed.")
    data <- data[!is.na(data[[idvar]]),]
  }
  
  # Number of clusters:
  nCluster <- length(unique(data[[idvar]]))
  
  # If beepvar is missing, add beepvar:
  if (missing(beepvar)){
    # First reorder by id:
    data <- data[order(data[[idvar]]),]
    
    beepvar <- 'BEEPVAR'
    data[[beepvar]] <- unlist(tapply(data[[idvar]],data[[idvar]],seq_along))
  } 
  # 
  # # Add dayvar:
  # if (missing(dayvar)){
  #   dayvar <- 'DAYVAR'
  #   data[[dayvar]] <- 1
  # } 
  # 
  # # Combine dayvar and beepvar:
  # alldays <- sort(unique(data[[dayvar]]))
  # allsubjects <- sort(unique(data[[idvar]]))
  # 
  # # Make first day first beeps:
  # scaleback <- function(x, start = 1) {
  #   x - min(x) + start
  # }
  # 
  # # For every subject fix the beeps:
  # for (p in allsubjects){
  #   # First day:
  #   ind1 <- data[[dayvar]] == alldays[1] & data[[idvar]] == allsubjects[p]
  #   if (any(ind1)){
  #     data[[beepvar]][ind1] <- scaleback(data[[beepvar]][ind1])  
  #   }
  #   
  #   
  #   # All other days:
  #   if (length(alldays)>1){
  #     for (i in 2:length(alldays)){
  #       
  #       ind_prev <- data[[dayvar]] %in% alldays[1:(i-1)] & data[[idvar]] == allsubjects[p]
  #       ind_cur <- data[[dayvar]] == alldays[i] & data[[idvar]] == allsubjects[p]
  #       if (any(ind_cur)){
  #         data[[beepvar]][ind_cur] <- scaleback(data[[beepvar]][ind_cur],max(data[[beepvar]][ind_prev])+1+nightskip)  
  #       }
  #       
  #       
  #     }
  #   }
  #   
  # }
  # browser()
  # else {
  # Check for integers:
  if (any(data[[beepvar]][na.omit(data[[beepvar]])] %% 1!=0)){
    stop("'beepvar' does not encode integer values.")
  }
  
  # Check for doubles:
  if (any(tapply(data[[beepvar]],data[[idvar]],function(x)any(duplicated(x)), simplify = TRUE))){
    stop("'beepvar' contains double values for one or more cases.")
  }
  
  # Overwrite minimum to always be 1:
  # data[[beepvar]] <- data[[beepvar]] - min(data[[beepvar]],na.rm = TRUE) + 1
  # }
  
  
  # Add group:
  if (missing(groups)){
    groups <- "GROUPID"
    data[[groups]] <- "fullsample"
  }
  
  # Max in cluster:
  maxInCluster <- max(data[[beepvar]])
  
  # Standardize the data:
  if (standardize == "z"){
    for (v in seq_along(vars)){
      data[,vars[v]] <- (data[,vars[v]] - mean(data[,vars[v]],na.rm=TRUE)) / sd(data[,vars[v]],na.rm=TRUE)
      # data[[vars[v]]] <- (data[[vars[v]]] - mean(data[[vars[v]]],na.rm=TRUE)) / sd(data[[vars[v]]],na.rm=TRUE)
    }
  } else if (standardize == "quantile"){
    for (v in seq_along(vars)){
      data[,vars[v]] <- quantiletransform(data[,vars[v]])
    }
  }
  
  
  # To long format:
  datalong <- tidyr::gather(data,variable,value,vars)
  
  # Transform data to wide format:
  datawide <- tidyr::pivot_wider(datalong, id_cols = c(idvar,groups),  values_from = "value", names_from = c("variable",beepvar))
  
  # Now make a design matrix:
  rowVars <- vars
  colVars <- as.character(seq(maxInCluster))
  design <- matrix(NA, length(rowVars), length(colVars))
  for (i in seq_along(rowVars)){
    for (j in seq_along(colVars)){
      varName <- paste0("^",rowVars[i],"_",colVars[j],"$")
      whichVar <- which(grepl(varName, names(datawide)))
      if (length(whichVar) == 1){
        design[i,j] <- paste0(rowVars[i],"_",colVars[j])
      }
    }
  }
  
  
  # Form model:
  model <- dlvm1(datawide, 
                 vars = design, 
                 groups = groups,
                 estimator = estimator,
                 ...
  )
  
  # Return:
  return(model)
}