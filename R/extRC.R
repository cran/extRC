extRC <-function(N, mod, k=min(dim(N))-1, la=0, marg.cons = c("free","equal","shift")){

# preliminaries
  marg.cons = match.arg(marg.cons)
  lev = dim(N); I = lev[1]; J = lev[2]
  y = c(t(N))
  out = MatIn(lev,mod)
  R0 = out$R0; R1 = out$R1; C0 = out$C0; C1 = out$C1
  J00 = out$J00; J01 = out$J01; J10 = out$J10; J11 = out$J11
  G = diag(I*J)[,-1]
  if(marg.cons=="free"){
    Cmg = matrix(0,0,I+J-2)
  }else{
    if(I!=J) stop("With constraints on the marginal probabilities, the number of columns must be equal to the number of rows")
    if(marg.cons=="equal") Cmg = cbind(diag(I-1),-diag(J-1))
    if(marg.cons=="shift") Cmg = dfm(I-1)%*%cbind(diag(I-1),-diag(J-1))
  }
  Ur = (dfm(I-1)%*%diag(c(1/(dfm(I)%*%(1:I)))))%x%cbind(1,matrix(0,1,J-2))
  Uc = cbind(1,matrix(0,1,I-2))%x%(dfm(J-1)%*%diag(c(1/(dfm(J)%*%(1:J)))))
  Cjn = matrix(0,0,(I-1)*(J-1))

# fit model
  np = length(la)
  if(np==1){
    Model = list(G=G,R0=R0,R1=R1,C0=C0,C1=C1,J00=J00,J01=J01,J10=J10,J11=J11,lev=lev,k=k,Cmg=Cmg,Cjn=Cjn,la=la)
    out = MainRC(y,Model,output=TRUE)
    P = matrix(out$pj,I,J,byrow=TRUE)
    dimnames(P) = list(X=1:I,Y=1:J)
    out$P = P
    etaX = out$eta[1:(I-1)]; etaY = out$eta[(I-1)+(1:(J-1))]
    Eta = matrix(out$eta[(I+J-2)+(1:(I-1)*(J-1))],I-1,J-1)
    dimnames(Eta) = list(X=1:(I-1),Y=1:(J-1))
    out$etaX = etaX; out$etaY = etaY; out$Eta = Eta
  }else{
    dev = rep(0,np)
    for(j in 1:np){
      Model = list(G=G,R0=R0,R1=R1,C0=C0,C1=C1,J00=J00,J01=J01,J10=J10,J11=J11,lev=lev,k=k,Cmg=Cmg,Cjn=Cjn,la=la[j])
      out = MainRC(y,Model)
      dev[j] = out$dev
    }
    out = list(la=la,dev=dev)
  }
  out$call = match.call()
  class(out) = "extRC"

# output
  return(out)

}