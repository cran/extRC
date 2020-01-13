cuby <- function(g){

# compute parameters
  X = rbind(c(1,0,0,0), # lk(0)
            c(0,1,0,0), # d1(0)
            c(1,1/4,1/16,1/64), # lk(1/4)
            c(0,1/2,1/4,1/8)) # lk(1/2)
  be = solve(t(X)%*%X)%*%t(X)%*%g
  
# find max
  c = be[2]; b = 2*be[3]; a = 3*be[4]
  dis = b^2-4*a*c; 
  if(dis<0 & a>=0){
      s = 0.2
  }else if(dis<0 && a<0){
      s = 0.05
    }else{
      dis = sqrt(dis) 
      xn = (-b-dis)/(2*a); xp=(-b+dis)/(2*a)
      Hn = 2*be[3]+6*be[4]*xn
      Hp = 2*be[3]+6*be[4]*xp
      if(xn>0 & Hn<0){
        s=xn
      }else{
        if(xp>0 & Hp<0){
          s = xp
        }else{
          s=0.05
        }
      }
    }

# output
  s = s

}