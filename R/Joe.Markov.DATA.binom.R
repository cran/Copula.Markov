Joe.Markov.DATA.binom  = function(n,size,prob,alpha){
  
  U_t = rep(NA,n)
  U_t[1] = runif(1)
  d_Joe = function(V){
    
    ( (1-U)^alpha + (1-V)^alpha - (1-U)^alpha*(1-V)^alpha )^(1/alpha-1) * 
      ( (1-U)^(alpha-1) - (1-V)^alpha * (1-U)^(alpha-1) ) - W
    
  }
  
  for( i in c(2:n) ){
    U = U_t[i-1]
    W = runif(1)
    
    U_t[i] = uniroot(d_Joe,interval=c(0.0000001,0.99999))$root
    
  }
  
  return(qbinom(p = U_t, size = size, prob = prob))
  
}





