tril <- function(M){

# preliminaries
  nr = nrow(M)
  nc = ncol(M)
  if(nr!=nc) stop("non-square matrix")
# put zeros in suitable elements
  for(i in 1:(nr-1)) M[i,(i+1):nr] = 0
# output
  N = M

}