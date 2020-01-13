Drank <- function(ga,lev,k,der=FALSE){

# preliminaries
  r = lev[1]-1; c = lev[2]-1; fr = ga
  if(k==min(r,c)){
    out = list(fr=NULL)
    if(der) out$Dfr=matrix(0,0,r*c)
    return(out)
  }else{
    Dfr = diag(r*c)
  }

# cicle across ranks
  for(it in 1:k){
    ind = which(abs(fr)>10^(-15))
    if(is.null(ind)){
      print('lack of rank')
      return()
    }else{
      ind = ind[1]; rc = r*c; irc = 1:rc; E = diag(rc)
      ro = ceiling(ind/c); co = ind-c*(ro-1)
      iro = c*(ro-1)+(1:c); ico = c*((1:r)-1)+co
      irc = setdiff(irc,union(iro,ico))
      iro = setdiff(iro,ind); ico = setdiff(ico,ind)
      f = fr[irc]-(fr[ico]%x%fr[iro])/fr[ind]
      if(der){
        D = E[ico,,drop=FALSE]%x%fr[iro]+fr[ico]%x%E[iro,,drop=FALSE]
        D = E[irc,]-D/fr[ind]
        D = D+(fr[ico]%x%fr[iro])%o%E[ind,]/fr[ind]^2
        Dfr = D%*%Dfr
      }
      r = r-1; c = c-1; fr = f
    }
  }

# output
  out = list(fr=fr)
  if(der) out$Dfr=Dfr
  return(out)

}