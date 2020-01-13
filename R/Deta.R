Deta <-
function(the,Model,der=FALSE){

# preliminaries
  G = Model$G; R0 = Model$R0; R1 = Model$R1; C0 = Model$C0; C1 = Model$C1
  J00 = Model$J00; J01 = Model$J01; J10 = Model$J10; J11 = Model$J11
  la = Model$la
  p = exp(G%*%the); p=p/sum(p)
  
# marginal logits
  Lr = log(R1%*%p)-log(R0%*%p); Lc = log(C1%*%p)-log(C0%*%p)

# interactions
  if(la==0){
    f11 = J11%*%p; f10 = J10%*%p; f01 = J01%*%p; f00 = J00%*%p
    Int = log(f11)-log(f10)-log(f01)+log(f00)
  }else{
    d11 = 1/((R1%*%p)%x%(C1%*%p)); d10 = 1/((R1%*%p)%x%(C0%*%p))
    d01 = 1/((R0%*%p)%x%(C1%*%p)); d00 = 1/((R0%*%p)%x%(C0%*%p))
    f11 = (J11%*%p)*d11; f10 = (J10%*%p)*d10; f01 = (J01%*%p)*d01; f00 = (J00%*%p)*d00
    Int = (f11^la-f10^la-f01^la+f00^la)/la
  }

# parameters
  eta = c(Lr,Lc,Int)

# Derivatives
  if(der){
    Mg = rbind((1/(R1%*%p))*R1/c(R1%*%p)-R0/c(R0%*%p),
               (1/(C1%*%p))*C1-(1/(C0%*%p))*C0)
    if(la==0){
           Dj = J11/f11-J10/f10-J01/f01+J00/f00
    }else{
      g11 = f11^(la-1); g10 = f10^(la-1); g01 = f01^(la-1); g00 = f00^(la-1)
      D11 = d11*J11-((J11%*%p)*d11^2)%*%(R1%x%(C1%*%p)+(R1%*%p)%x%C1)
      D10 = d10*J10-((J10%*%p)*d10^2)%*%(R1%x%(C0%*%p)+(R1%*%p)%x%C0)
      D01 = d01*J01-((J01%*%p)*d01^2)%*%(R0%x%(C1%*%p)+(R0%*%p)%x%C1)
      D00 = d00*J00-((J00%*%p)*d00^2)%*%(R0%x%(C0%*%p)+(R0%*%p)%x%C0)
      Dj = g11*D11-g10*D10-g01*D01+g00*D00
    }
    Der = rbind(Mg,Dj)*(p*G-(p%o%p)%*%G);
  }

# output
    out = list(eta=eta)
    if(der) out$Der=Der
    return(out)

}
