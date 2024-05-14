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
    data <- data %>% dplyr::group_by(.data[[dayvar]],.data[[idvar]]) %>% 
      dplyr::mutate(BEEP = seq_len(n()))
  }
  
  # Vars:
  if (is.null(vars)){
    vars <- names(data[!names(data)%in%c(idvar,dayvar,beepvar)])
  }
  
  
  # Only retain important columns:
  data <- data[,c(vars,idvar,dayvar,beepvar)]
  
  # Center and scale data:
  for (v in vars){
    data[,v] <- as.numeric(scale(data[,v], center, scale))
  }
  

  # Obtain person specific means:
 MeansData <-data %>% dplyr::group_by(.data[[idvar]]) %>% dplyr::summarise_at(funs(mean(.,na.rm=TRUE)),.vars = vars)
 
  #FIXME: Replace with?
 data %>% dplyr::group_by(.data[[idvar]]) %>% dplyr::summarise_at(list(~mean(.,na.rm=TRUE)),.vars = vars)
 
  
  # Within-person center:
  if (centerWithin){
    # Only if N > 1 (very minimal floating point error can lead to different layout to older version otherwise)
    if (length(unique(data[[idvar]])) > 1){
      data <- data %>% dplyr::group_by(.data[[idvar]]) %>% dplyr::mutate_at(funs(scale(.,center=TRUE,scale=FALSE)),.vars = vars)         
    }
  }

  # From mlVAR: Augment data:
  # Augment the data
  augData <- data
  
 
  
  # Check for errors in data:
  beepsummary <- data %>% group_by(.data[[idvar]],.data[[dayvar]],.data[[beepvar]]) %>% tally
  if (any(beepsummary$n!=1)){
    print_and_capture <- function(x)
    {
      paste(capture.output(print(x)), collapse = "\n")
    }
    
    warning(paste0("Some beeps are recorded more than once! Results are likely unreliable.\n\n",print_and_capture(
      beepsummary %>% filter(.data[["n"]]!=1) %>% select(.data[[idvar]],.data[[dayvar]],.data[[beepvar]]) %>% as.data.frame
    )))
  }
  
  beepsPerDay <-  dplyr::summarize(data %>% group_by(.data[[idvar]],.data[[dayvar]]), 
                                                    first = min(.data[[beepvar]],na.rm=TRUE),
                                                    last = max(.data[[beepvar]],na.rm=TRUE))
 
  # all beeps:
  allBeeps <- expand.grid(unique(data[[idvar]]),unique(data[[dayvar]]),seq(min(data[[beepvar]],na.rm=TRUE),max(data[[beepvar]],na.rm=TRUE))) 
  names(allBeeps) <- c(idvar,dayvar,beepvar)
  
  # Left join the beeps per day:
  allBeeps <- allBeeps %>% dplyr::left_join(beepsPerDay, by = c(idvar,dayvar)) %>% 
      dplyr::group_by(.data[[idvar]],.data[[dayvar]]) %>% dplyr::filter(.data[[beepvar]] >= .data$first, .data[[beepvar]] <= .data$last)%>%
      dplyr::arrange(.data[[idvar]],.data[[dayvar]],.data[[beepvar]])
  
  
  # Enter NA's:
  augData <- augData %>% dplyr::right_join(allBeeps, by = c(idvar,dayvar,beepvar)) %>%
    arrange(.data[[idvar]],.data[[dayvar]],.data[[beepvar]])
  
  # Obtain data_c (slice away first row per day/subject):
  data_c <- augData %>% ungroup %>% dplyr::select(all_of(vars))
  
  # Lagged datasets:
  data_l <- do.call(cbind,lapply(lags, function(l){
    data_lagged <- augData %>% dplyr::group_by(.data[[idvar]],.data[[dayvar]]) %>% dplyr::mutate_at(funs(shift),.vars = vars) %>% ungroup %>% dplyr::select(all_of(vars))
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