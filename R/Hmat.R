Hmat <-function(G){

# matrix transformation
  t = nrow(G); sg = colSums(G)
  u = rep(1,t)
  H = solve(t(G)%*%G-sg%o%sg/t)%*%(t(G)-sg%o%u/t)

# output
  H = H

}