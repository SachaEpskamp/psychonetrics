# this function will automatically identify models:
identify <- function(x){
  if (x@model == "ggm" | x@model == "precision"| x@model == "gvar" | x@model == "varcov" | x@model == "cholesky"){
    # Nothing to do..
    return(x)
  }

  if (x@model == "lnm"){
    x <- identify_lnm(x)
    return(x)
  } else if (x@model == "rnm"){
    x <- identify_rnm(x)
    return(x)
  } else if (x@model == "lvm"){
    x <- identify_lvm(x)
    return(x)
  }
  
  stop("Model not supported...")
}
