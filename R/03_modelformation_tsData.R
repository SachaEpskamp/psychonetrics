# Mostly copied from graphicalVAR...

shift <- function(x,n){
  length <- length(x)
  c(rep(NA,n),x)[1:length]
}

# Data preperation function:
tsData <- function(data,
                     vars =NULL,
                     beepvar = NULL,
                     dayvar = NULL,
                     idvar = NULL,
                   groupvar = NULL,
                   lags = 1,
                     scale = FALSE,
                   center = FALSE,
                     centerWithin = FALSE # False if idvar is missing, true otherwise
                   ){
  # If not missing groupvar, just do this per group:

  if (!is.null(groupvar)){
    groups <- unique(data[[groupvar]])
    res <- lapply(seq_along(groups),function(g){
      dat <- data[data[[groupvar]] == groups[g],names(data)!=groupvar] 
      dat <- tsData(dat, vars, beepvar, dayvar, idvar)
      dat[[groupvar]] <- groups[g]
      dat
    })
    return(do.call(rbind,res))
  }
  
  
  . <- NULL
  deleteMissings = FALSE
  
  data <- as.data.frame(data)
  
  # defaults:
  # scale <- FALSE
  # centerWithin <- TRUE
  # 
  # Add subject:
  if (is.null(idvar)){
    idvar <- "ID"
    data[[idvar]] <- 1
  }
  
  # Add day:
  if (is.null(dayvar)){
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  }

  # Add beepvar:
  if (is.null(beepvar)){
    beepvar <- "BEEP"
    data <- data %>% dplyr::group_by_(dayvar,idvar) %>% 
      dplyr::mutate_(BEEP = ~seq_len(n()))
  }
  
  # Vars:
  if (is.null(vars)){
    vars <- names(data[!names(data)%in%c(idvar,dayvar,beepvar)])
  }
  
  
  # Only retain important columns:
  data <- data[,c(vars,idvar,dayvar,beepvar)]
  
  # Center and scale data:
  data[,vars] <- scale(data[,vars], center, scale)
  
  # Obtain person specific means:
  MeansData <- data %>% dplyr::group_by_(idvar) %>% dplyr::summarise_at(funs(mean(.,na.rm=TRUE)),.vars = vars)
  
  # Within-person center:
  if (centerWithin){
    # Only if N > 1 (very minimal floating point error can lead to different layout to older version otherwise)
    if (length(unique(data[[idvar]])) > 1){
      data <- data %>% dplyr::group_by_(idvar) %>% dplyr::mutate_at(funs(scale(.,center=TRUE,scale=FALSE)),.vars = vars)          
    }
  }

  # From mlVAR: Augment data:
  # Augment the data
  augData <- data
  
  # Add missing rows for missing beeps
  beepsPerDay <-  eval(substitute(dplyr::summarize_(data %>% group_by_(idvar,dayvar), 
                                                    first = ~ min(beepvar,na.rm=TRUE),
                                                    last = ~ max(beepvar,na.rm=TRUE)), 
                                  list(beepvar = as.name(beepvar))))
  
  # all beeps (with one extra beep before each measurement of the first day:
  allBeeps <- expand.grid(unique(data[[idvar]]),unique(data[[dayvar]]),seq(min(data[[beepvar]],na.rm=TRUE)-1,max(data[[beepvar]],na.rm=TRUE))) 
  names(allBeeps) <- c(idvar,dayvar,beepvar)
  
  # Left join the beeps per day:
  allBeeps <- eval(substitute({
    allBeeps %>% dplyr::left_join(beepsPerDay, by = c(idvar,dayvar)) %>% 
      # dplyr::group_by_(idvar,dayvar) %>% dplyr::filter_(~BEEP >= first, ~BEEP <= last)%>%
      dplyr::group_by_(idvar,dayvar) %>% dplyr::filter_(~BEEP >= first-1, ~BEEP <= last)%>%
      dplyr::arrange_(idvar,dayvar,beepvar)
  },  list(BEEP = as.name(beepvar))))

  
  # Enter NA's:
  augData <- augData %>% dplyr::right_join(allBeeps, by = c(idvar,dayvar,beepvar))
  
  # Obtain data_c (slice away first row per day/subject):
  data_c <- augData %>% ungroup %>% dplyr::select_(.dots = vars)#  %>% dplyr::group_by_(idvar,dayvar) %>% dplyr::slice(-1)
  
  # Lagged datasets:
  data_l <- do.call(cbind,lapply(lags, function(l){
    data_lagged <- augData %>% dplyr::group_by_(idvar,dayvar) %>% dplyr::mutate_at(funs(shift),.vars = vars) %>% ungroup %>% dplyr::select_(.dots=vars)
    names(data_lagged) <- paste0(vars,"_lag",l)
    data_lagged
  }))
  
  # Remove missing (full row):
  isNA <- rowSums(is.na(data_c)) == ncol(data_c)
  data_c <- data_c[!isNA,]
  data_l <- data_l[!isNA,]
  
  # Combine them:
  fulldata <- as.data.frame(cbind(data_l,data_c))
  return(fulldata)
  
  # data_l <- augData %>% dplyr::group_by_(idvar,dayvar) %>% dplyr::slice(-n())
# 
#   
  # # Remove rows with missings:
#   if (deleteMissings){

#   }

  # Return datasets:
  # Results <- list(
  #   data = augData,
  #   data_c = data_c[,vars],
  #   data_l = cbind(1,data_l),
  #   data_means = MeansData,
  #   vars=vars,
  #   idvar=idvar,
  #   dayvar=dayvar,
  #   beepvar=beepvar,
  #   lags = lags
  # )
  # 
  # class(Results) <- "tsData"
  # return(Results)
}