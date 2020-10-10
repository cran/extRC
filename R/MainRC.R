MainRC <-function(y,Model,the0=NULL,output=FALSE){

# Preliminaries
  G = Model$G; la = Model$la
  eps = 10^(-8); tol = 10^(-8); maxit = 600
  n = sum(y); y = y+eps*(y==0); y = y*n/sum(y)

# Raw estimate  
  pj = y/n; LKu = c(y%*%log(pj)); H = Hmat(G)
  if(is.null(the0)) the0 = H%*%log(pj)/5
  out = PraD(the0,Model,der=TRUE)
  hdis = out$hdis; Hdis = out$Hdis

# Iterate
  it = 0; fc = rep(0,4); tmax = 200
  if(sum(abs(hdis)) <= tol){
    th1 = the0
    dev = 0
  }
  while(sum(abs(hdis)) > tol & it<maxit){
    it=it+1
    pj = c(exp(G%*%the0)); pj = pj/sum(pj)
# step direction
    X = Null(t(Hdis)); hi = ginv(Hdis)%*%hdis; sc = t(G)%*%(y-n*pj)
    FF = n*t(G)%*%(pj*G-(pj%o%pj)%*%G)
    v0 = hi+solve(FF)%*%sc; A = t(X)%*%FF%*%X
    a = t(X)%*%FF%*%v0; de = c(X%*%solve(A)%*%a-hi)
# try optimize step length (3 trials)
    fc[1] = y%*%log(pj)/n-hdis%*%hdis/2
    fc[2] = de%*%c(t(G)%*%(y/n-pj)-t(Hdis)%*%hdis)
    th = the0+de/4; pj = c(exp(G%*%th)); pj = pj/sum(pj)
    th = pmax(pmin(th,tmax),-tmax)
    out = PraD(th,Model,der=TRUE)
    hdis = out$hdis; Hdis = out$Hdis
    fc[3] = y%*%log(pj)/n-hdis%*%hdis/2; th = the0+de/2 
    pj = c(exp(G%*%th)); pj = pj/sum(pj); hdis = PraD(th,Model)$hdis
    fc[4] = y%*%log(pj)/n-hdis%*%hdis/2; st = cuby(fc) 
# reconstruct
    th1 = the0+st*de; pj = c(exp(G%*%th1)); pj = pj/sum(pj)
    th1 = pmax(pmin(th1,tmax),-tmax)
    out = PraD(th1,Model)
    eta = out$eta; hdis = out$hdis
    LK = c(y%*%log(pj)); the0 = th1
# Computes deviance
    dev = 2*(LKu-LK)
  }
  df = length(hdis)
# output
  out = list(dev=dev,df=df,pj=pj,it=it,dis=sum(abs(hdis)))
  if(output) out$eta = eta 
  return(out)

}