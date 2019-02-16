# Make a list with for each matrix two elements: est containing the current estimate, and par containing parameter numbers (0 = fixed)
addMatrix <- function(x,nPars = 0, symmetrical = FALSE,start=0){
  matList <- list(
    cur = ifelse(is.na(x),start,x),
    par = ifelse(is.na(x),TRUE,FALSE)
  )
  if (symmetrical){
    matList$par[matList$par & lower.tri(matList$par,diag=TRUE)] <- nPars + seq_len(sum(matList$par & lower.tri(matList$par,diag=TRUE))) # Add par numbers
    matList$par[matList$par & upper.tri(matList$par)] <- t(matList$par)[matList$par & upper.tri(matList$par)]
  } else {
    matList$par[matList$par] <- nPars + seq_len(sum(matList$par)) # Add par numbers
  }
  matList$par <- 1 * matList$par # Make numeric in case there was no par number
  matList
}

# Helper function to delete unnesisary columns:
deleteCols <- function(x, mat, parList, symmetrical=FALSE){
  if (!symmetrical){
    x <- x[,parList[[mat]]$par != 0]
  } else {
    x <- x[,parList[[mat]]$par[lower.tri(parList[[mat]]$par,diag=TRUE)] != 0]
  }
  x
}

# Function to update the parList (used later)
updateParList <- function(x,parList){
  for (i in seq_along(parList)){
    if (sum(parList[[i]]$par!=0)>0){
      parList[[i]]$cur[parList[[i]]$par!=0] <- x[parList[[i]]$par[parList[[i]]$par!=0]]
    }
  }
  parList
}

# Funtion to form matrices:
formMatrix <- function(x,mat,parList){
  res <- parList[[mat]]$cur
  if (max(parList[[mat]]$par)!=0){
    res[parList[[mat]]$par!=0] <- x[parList[[mat]]$par[parList[[mat]]$par!=0]]    
  }
  res
}
