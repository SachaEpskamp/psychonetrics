# Inner function:
replicator_inner_multicore <- function(x,input,expr,packages = c("psychonetrics","magrittr")){
  lapply(packages, library, character.only=TRUE)
  suppressWarnings(eval(expr))
}
replicator_inner_singlecore <- function(x,input,expr){
  suppressWarnings(eval(expr))
}

# Function to replicate a call:
replicator <- function(
  input, # some input
  expression, # some expression
  results = c("parameters","full","matrices"), # returns to full if results are not psychonetrics
  reps = 10,
  nCores = 1,
  verbose = TRUE,
  env = parent.frame(n=1),
  packages = c("psychonetrics","magrittr"), # Packages to load
  export
){
  results <- match.arg(results)
  
  # Substitute expression:
  # FIXME: Ugly way...
  expr <- parse(text=paste0("input %>% ",deparse(substitute(expression))))
  
  
  # Ncores > 1:
  if (nCores > 1){
    nClust <- nCores - 1
    cl <- makePSOCKcluster(nClust)  
    if (!missing(export)){
      parallel::clusterExport(cl, export)  
    }
  # Run:   
    #FIXME: At the moment I copy the global workspace, BAD SOLUTION
    Results <- pblapply(seq_len(reps),FUN = replicator_inner_multicore, input=input, expr=expr, packages=packages,cl = cl)
    
    # Stop the cluster:
    stopCluster(cl)
  } else {
    Results <- pblapply(seq_len(reps),FUN = replicator_inner_singlecore, input=input,  expr=expr)
  }

  if (results == "full"){
    return(Results)
  } else if (results == "parameters") {
    if (!is(Results[[1]], "psychonetrics")){
      warning("Replications did not return psychonetrics objects, returning raw replications instead.")
      return(Results)
    } 
    # Obtain parameter table:
    pars <- do.call(rbind,lapply(Results,function(x)x@parameters$est))
    pars <- as.data.frame(pars)
    names(pars) <- paste0(Results[[1]]@parameters$matrix,"[",Results[[1]]@parameters$row,",",Results[[1]]@parameters$col,"]")
    return(pars)
  } else if (results == "matrices"){
    if (!is(Results[[1]], "psychonetrics")){
      warning("Replications did not return psychonetrics objects, returning raw replications instead.")
      return(Results)
    } 
    
    
    Matrices <- lapply(Results, function(x) x@modelmatrices)
    matNames <- names(Matrices[[1]][[1]])
    nGroups <- length(Matrices[[1]])
    outList <- vector("list",nGroups)
    for (g in 1:nGroups){
      for (mat in matNames){
        outList[[g]][[mat]] <- lapply(Matrices,function(x)x[[g]][[mat]])
      }
    }
    return(outList)
  }
  

}