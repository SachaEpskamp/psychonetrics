goodNum <- function(x){
  sapply(x,function(xx){
    if (is.na(xx))return("")
    if (xx < 0.0001){
      return("< 0.0001")
    }
    digits <- max(0,floor(log10(abs(xx))) + 1)
    isInt <- xx%%1 == 0
    gsub("\\.$","",formatC(signif(unlist(xx),digits=digits+(!isInt)*2), digits=digits+(!isInt)*2,format="fg", flag="#"))
  })  
}



goodNum2 <- function(x){
  sapply(x,function(xx){
    if (is.na(xx))return("")
    if (xx < 0.0001 & xx > -0.0001){
      return("~0")
    }
    digits <- max(0,floor(log10(abs(xx))) + 1)
    isInt <- xx%%1 == 0
    gsub("\\.$","",formatC(signif(unlist(xx),digits=digits+(!isInt)*2), digits=digits+(!isInt)*2,format="fg", flag="#"))
  })  
}