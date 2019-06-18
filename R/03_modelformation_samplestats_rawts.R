# # Inner function to make sample stats object:
# samplestats_rawts <- function(
#   data, # Dataset
#   vars, # character indicating the variables Extracted if missing from data - group variable
#   groups, # ignored if missing. Can be character indicating groupvar, or vector with names of groups
#   ... # Ignored arguments
# ){
#   # Check data:
#   if (missing(data) ){
#     stop("'data' may not be missing when rawts = TRUE")
#   }
#   
#   # If data is supplied:
#   if (!is.data.frame(data) & !is.matrix(data)){
#     stop("'data' must be a data frame or matrix")
#   }
#   if (is.matrix(data)){
#     data <- as.data.frame(data)
#   }
#   # If group variable is missing, add (dummy):
#   if (missing(groups)|| is.null(groups)){
#     groups <- "singlegroup"
#     data[[groups]] <- "singlegroup"
#   }
#   # Extract group names:
#   groupNames <- unique(data[[groups]])
#   
#   # number of groups:
#   nGroup <- length(groupNames)
#   
#   # Overwrite group with integer:
#   data[[groups]] <- match(data[[groups]], groupNames)
#   
#   # If vars is missing, obtain from data:
#   if (missing(vars)){
#     vars <- names(data[,names(data)!=groups])
#   }
# 
#   # Number of variables:
#   nVars <- length(vars)
#   
#   
#   # Create covs and means arguments:
#   if (nGroup == 1){
#     vecData <- as.vector(t(as.matrix(data[,vars])))
#     vecData_noNA <- na.omit(vecData) # The missingness pattern is handled later
#     covs <- cors <- list(Matrix(0,nrow=length(vecData_noNA),ncol=length(vecData_noNA)))
#     means <- list(Matrix(vecData_noNA))
#     groupNames <- unique(data[[groups]])
#     missingness <- list(Matrix(is.na(data[,vars])))
#     
#   } else {
#     covs <- list()
#     cors <- list()
#     means <- list()
#     missingness <- list()
#     groupNames <- unique(data[[groups]])
#     
#     for (g in 1:nGroup){
#       subData <- data[data[[groups]] == g,c(vars)]
#       vecData <- as.vector(t(as.matrix(subData)))
#       vecData_noNA <- na.omit(vecData) # The missingness pattern is handled later
#       covs[[g]] <- cors[[g]] <- Matrix(0,nrow=length(vecData_noNA),ncol=length(vecData_noNA))
#       means[[g]] <- Matrix(vecData_noNA)
#       missingness[[g]] <- Matrix(is.na(data[,vars]))
#     }
#     
#   }
#   # nobs <- as.vector(tapply(data[[groups]],data[[groups]],length))
#   
#   # Set names:
#   names(covs) <- groupNames
#   names(means) <- groupNames
#   
#   # Generate samplestats object:
#   object <- generate_psychonetrics_samplestats(covs = covs, cors = cors, means = means)
#   
#   # Fill groups:
#   object@groups <- data.frame(
#     label = groupNames,
#     id = seq_along(groupNames),
#     nobs = 1, stringsAsFactors = FALSE
#   )
#   
#   # Fill variables:
#   object@variables <- data.frame(
#     label = vars,
#     id = seq_along(vars), stringsAsFactors = FALSE
#   )
#   
#   # Fill missingness:
#   object@missingness <- missingness
# 
#   
#   # Return object:
#   return(object)
# }