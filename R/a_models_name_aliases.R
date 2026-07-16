# Primary short names (0.16.7 naming cleanup): the "1" suffix in several model
# names referred to the lag order (VAR(1)), but no lag-2 frameworks are planned,
# so the plain names are now the primary documented names. The *1 names remain
# fully supported aliases (no deprecation). var1()/gvar() are unchanged:
# exporting var() would mask stats::var, and "var1" reads naturally as VAR(1).

dlvm <- function(...){
  dlvm1(...)
}

tsdlvm <- function(...){
  tsdlvm1(...)
}

meta_var <- function(...){
  meta_var1(...)
}
