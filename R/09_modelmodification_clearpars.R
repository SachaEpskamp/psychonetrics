# Clears parameter indices (input is the parameter table)
clearpars <- function(x, ind){

  x$std[ind] <- NA
  x$se[ind] <- NA
  x$p[ind] <- NA
  x$se_boot[ind] <- NA
  x$p_boot[ind] <- NA
  x$mi[ind] <- NA
  x$pmi[ind] <- NA
  x$epc[ind] <- NA
  x$mi_free[ind] <- NA
  x$pmi_free[ind] <- NA
  x$epc_free[ind] <- NA
  x$mi_equal[ind] <- NA
  x$pmi_equal[ind] <- NA
  x$epc_equal[ind] <- NA
  
  
  x
}