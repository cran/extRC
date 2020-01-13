plot.extRC <-function(x, ...){

# preliminaries
  out = x
# plot output
  plot(out$la,out$dev,type="l",xlab="lambda",ylab="deviance")

}
