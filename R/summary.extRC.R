summary.extRC <- function(object, ...){

# preliminaries
  out = object
# print output
  cat("\nCall:\n")
  print(out$call)
  if(length(out$dev)==1){
    cat("\ndeviance:\n")
    print(round(out$dev,2))
    cat("\ndf:\n")
    print(out$df)
    cat("\njoint probability matrix:\n")
    print(round(out$P,4))
    cat("\nrow marginal logits:\n")
    print(round(out$etaX,3))
    cat("\ncolumn marginal logits:\n")
    print(round(out$etaY,3))
    cat("\nlog-odds ratios:\n")
    print(round(out$Eta,3))
    cat("\niteration number:\n")
    print(out$it)
    cat("\ndiscrepancy:\n")
    print(out$dis)
    cat("\n")
  }

}
