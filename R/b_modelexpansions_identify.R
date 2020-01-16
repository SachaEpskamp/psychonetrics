# this function will automatically identify models:
identify <- function(x){
  if (x@model == "ggm" | x@model == "precision"| x@model == "gvar" | x@model == "varcov" | x@model == "cholesky" | x@model == "var1" | x@model == "panelvar1" | x@model == "meta_varcov"){
    # Nothing to do..
    return(x)
  }

  # if (x@model == "lnm"){
  #   x <- identify_lnm(x)
  #   return(x)
  # } else if (x@model == "rnm"){
  #   x <- identify_rnm(x)
  #   return(x)
  # } else 
  if (x@model == "lvm"){
    x <- identify_lvm(x)
    return(x)
  } else if (x@model == "dlvm1"){
    x <- identify_dlvm1(x)
    return(x)
  } else if (x@model == "tsdlvm1"){
    x <- identify_tsdlvm1(x)
    return(x)
  } else if (x@model == "Ising"){
    x <- identify_Ising(x)
    return(x)
  }
  
  
  stop("Model not supported...")
}
