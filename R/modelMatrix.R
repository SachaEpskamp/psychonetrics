# Form the model matrix M :
Mmatrix <- function(x){
  M <- sparseMatrix(i = seq_len(nrow(x))[x$par!=0], j = x$par[x$par!=0], dims = c(nrow(x),max(x$par)))
  M
}