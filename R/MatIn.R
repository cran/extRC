MatIn <- function(lev,mod){

# Preliminaries
  Rd = tril(matrix(1,lev[1],lev[1]))
  Ru = t(Rd)
  Cd = tril(matrix(1,lev[2],lev[2]))
  Cu = t(Cd)
  Id = diag(max(lev))
  rd = 1:lev[1]-1; ru = 2:lev[1]
  cd = 1:lev[2]-1; cu = 2:lev[2] 
  
# Rows
  if(mod[1] == "g"){ # global
    R0 = Rd[rd,]; R1 = Ru[ru,]
  }
  if(mod[1] == 'c'){ # continuation
    R0=Id[rd,1:lev[1]]; R1 = Ru[ru,]
  }
  if(mod[1] == 'l'){ # adjacent
    R0=Id[rd,1:lev[1]]; R1 = Id[ru,1:lev[1]]     
  }

  # Cols
  if(mod[2] == "g"){ # global
    C0 = Cd[cd,]; C1 = Cu[cu,]
  }
  if(mod[2] == "c"){ # continuation
    C0 = Id[cd,1:lev[2]]; C1 = Cu[cu,]
  }
  if(mod[2] == "l"){ # adjacent
    C0 = Id[cd,1:lev[2]]; C1 = Id[cu,1:lev[2]]
  }

# Joint
  J00 = R0%x%C0; J01 = R0%x%C1
  J10 = R1%x%C0; J11 = R1%x%C1

# marginals
  R0 = R0%x%matrix(1,1,lev[2]); R1 = R1%x%matrix(1,1,lev[2])
  C0 = matrix(1,1,lev[1])%x%C0; C1 = matrix(1,1,lev[1])%x%C1

# output
  out = list(R0=R0,R1=R1,C0=C0,C1=C1,J00=J00,J01=J01,J10=J10,J11=J11)
  
}