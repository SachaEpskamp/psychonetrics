# Function to aggregate bootstraps:
aggregate_bootstraps <- function(
    sample,
    bootstraps
){
  # If not missing sample, check if it is a psychonetrics object:
  if (!missing(sample)){
    stopifnot(is(sample,"psychonetrics"))
  }
  
  # # obtain dots:
  # dots <- c(list(...),bootstraps)
  # 
  # # Check if any are lists:
  # is_list <- sapply(dots, function(x) is.list(x))
  # 
  # # Expand:
  # new_dots <- c(dots[!is_list], unlist(dots[is_list], recursive = FALSE))
  new_dots <- bootstraps
  
  # Now check if all are psychonetrics bootstap objects:
  all_good <- all(sapply(new_dots,function(x){
    if (!is(x,"psychonetrics")){
      return(FALSE)
    } else 
      return(x@sample@bootstrap)
  }))
  
  # check:
  if (!all_good){
    stop("Input 'bootstraps' is not a list of 'psychonetrics' objects that are bootstrapped.")
  }
  
  # Next check if the size of the model is the same for each bootstrap:
  n_pars <- sapply(new_dots,function(x) nrow(x@parameters))
  models <- sapply(new_dots,function(x) x@model)
  submodels <- sapply(new_dots,function(x) x@submodel)
  distributions <- sapply(new_dots,function(x) x@distribution)
  
  # Check if all equal:
  equal_test <- function(x){
    length(unique(x)) == 1
  }
  if (!equal_test(n_pars) || !equal_test(models) || !equal_test(submodels)  || !equal_test(distributions)){
    stop("Not all bootstraps are based on the same model.")
  }
  
  # Check if bootstrap parameters are the same:
  boot_sub <- sapply(new_dots,function(x) x@sample@boot_sub)
  boot_resample <- sapply(new_dots,function(x) x@sample@boot_resample)
  
  if (!equal_test(boot_sub) || !equal_test(boot_resample)){
    stop("Bootstrap parameters ('boot_sub' and 'boot_resample') are not equal across bootstrap samples.")
  }
  
  # Empty aggregate model:
  agg_mod <- generate_psychonetrics_bootstrap()
  
  # Overwrite info from first model:
  agg_mod@model <- models[1]
  agg_mod@submodel <- submodels[1]
  agg_mod@distribution <- distributions[1]
  agg_mod@verbose <- new_dots[[1]]@verbose
  
  agg_mod@boot_sub <- boot_sub[1]
  agg_mod@boot_resample <- boot_resample[1]
  
  agg_mod@matrices <- new_dots[[1]]@matrices
  agg_mod@types <- new_dots[[1]]@types
  # Store all models:
  agg_mod@models <- new_dots
  
  ## Form the parameter table:
  pars <- lapply(new_dots,function(x) x@parameters)
  
  # Add bootstrap:
  for (i in seq_along(pars)){
    pars[[i]]$boot <- i
    pars[[i]]$parrow <- seq_len(nrow(pars[[i]]))
  }
  
  # Combine:
  pars <- dplyr::bind_rows(pars)
  
  minNA <- function(x,na.rm=TRUE){
    if (length(x)==0) return(NA) else return(min(x,na.rm=na.rm))
  }
  
  
  maxNA <- function(x,na.rm=TRUE){
    if (length(x)==0) return(NA) else return(max(x,na.rm=na.rm))
  }
  
  # Aggregate:
  boot_pars <- pars %>% group_by(.data[["parrow"]],.data[["var1"]],.data[["var1_id"]],.data[["op"]],.data[["var2"]],.data[["var2_id"]],.data[["matrix"]],.data[["row"]],.data[["col"]],.data[["group"]],.data[["group_id"]]) %>% 
    summarize(
      avg = mean(.data[["est"]], na.rm=TRUE),
      se_boot = sd(.data[["est"]], na.rm=TRUE),
      min = minNA(.data[["est"]], na.rm=TRUE),
      q1 = quantile(.data[["est"]], 1/100, na.rm=TRUE),
      q2.5 = quantile(.data[["est"]], 2.5/100, na.rm=TRUE),
      q5 = quantile(.data[["est"]], 5/100, na.rm=TRUE),
      median = median(.data[["est"]], na.rm=TRUE),
      q95 = quantile(.data[["est"]], 95/100, na.rm=TRUE),
      q97.5 = quantile(.data[["est"]], 97.5/100, na.rm=TRUE),
      q99 = quantile(.data[["est"]], 99/100, na.rm=TRUE),
      max = maxNA(.data[["est"]], na.rm=TRUE),
      
      avg_non0 = mean(.data[["est"]][.data[["est"]]!=0], na.rm=TRUE),
      se_boot_non0 = sd(.data[["est"]][.data[["est"]]!=0], na.rm=TRUE),
      prop_non0 = mean(.data[["est"]]!=0, na.rm=TRUE),
      prop_non0_pos = mean(.data[["est"]]>0, na.rm=TRUE),
      prop_non0_neg = mean(.data[["est"]]<0, na.rm=TRUE),
      min_non0 = minNA(.data[["est"]][.data[["est"]]!=0], na.rm=TRUE),
      q1_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 1/100, na.rm=TRUE),
      q2.5_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 2.5/100, na.rm=TRUE),
      q5_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 5/100, na.rm=TRUE),
      median_non0 = median(.data[["est"]][.data[["est"]]!=0], na.rm=TRUE),
      q95_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 95/100, na.rm=TRUE),
      q97.5_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 97.5/100, na.rm=TRUE),
      q99_non0 = quantile(.data[["est"]][.data[["est"]]!=0], 99/100, na.rm=TRUE),
      max_non0 = maxNA(.data[["est"]][.data[["est"]]!=0], na.rm=TRUE)
    ) %>% ungroup %>%  arrange(.data[["parrow"]]) %>% select(-.data[["parrow"]])
  
  # Add sample if not missing:
  if (!missing(sample)){
    boot_pars$est_sample <- sample@parameters$est
    boot_pars$fixed_sample <- sample@parameters$fixed
  } else {
    boot_pars$est_sample <- NA
    boot_pars$fixed_sample <- NA
  }
  
 # Store in parameters:
  agg_mod@parameters <- bind_rows(agg_mod@parameters, boot_pars)
  

  
  # Aggregate fit measures:
  fit_measures <- lapply(new_dots,function(x) data.frame(
    measure = names(x@fitmeasures),
    value = unlist(x@fitmeasures))
  )
  for (i in seq_along(fit_measures)){
    fit_measures[[i]]$boot <- i
  }
  fit_measures <- bind_rows(fit_measures)
  
  # Aggregate:
  rownames(fit_measures) <- NULL
  fit_measures <- fit_measures %>% group_by(.data[["measure"]]) %>% 
    summarize(
     avg = mean(.data[["value"]], na.rm=TRUE),
     sd = sd(.data[["value"]], na.rm=TRUE),
     min = min(.data[["value"]], na.rm=TRUE),
     q1 = quantile(.data[["value"]], 1/100, na.rm=TRUE),
     q2.5 = quantile(.data[["value"]], 2.5/100, na.rm=TRUE),
     q5 = quantile(.data[["value"]], 5/100, na.rm=TRUE),
     median = median(.data[["value"]], na.rm=TRUE),
     q90 = quantile(.data[["value"]], 90/100, na.rm=TRUE),
     q95 = quantile(.data[["value"]], 95/100, na.rm=TRUE),
     q99 = quantile(.data[["value"]], 99/100, na.rm=TRUE),
     max = max(.data[["value"]], na.rm=TRUE)
    ) 
  
  agg_mod@fitmeasures <- fit_measures
  
  # Return:
  return(agg_mod)
  
}
