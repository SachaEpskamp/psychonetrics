emergencystart <- function(x){
  # Try fixing:
  # For all matrices in the model:
  allMats <- x@matrices[!grepl("(mu)|(tau)",x@matrices$name),]
  for (mat in seq_len(nrow(allMats))){
    # Adjust starting values. Identity matrix for posdef matrices, Near-null matrix for anything else:
    if (allMats$posdef[mat]){
      x@parameters$est[x@parameters$matrix == allMats$name[mat]] <- 
        ifelse(x@parameters$row[x@parameters$matrix == allMats$name[mat]] ==
                 x@parameters$col[x@parameters$matrix == allMats$name[mat]], 1, 0)
    } else {
      x@parameters$est[x@parameters$matrix == allMats$name[mat]] <- 
        (x@parameters$est[x@parameters$matrix == allMats$name[mat]] != 0) * 1e-7 * sign(x@parameters$est[x@parameters$matrix == allMats$name[mat]])
    }
    
  }
  
  return(x)
}