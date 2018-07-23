Clayton.Markov.DATA.binom <-
function(n,size,prob,alpha){

  ###############Clayton copula(to the pwer of -1/a) 
  C0=function(u1,u2,a)
  {
    {x=((1/(u1^(a))+1/(u2^(a))-1)) }
    x[x<=0]=10^(-200)
    z=x^(-1/a)
    z
  }
  
  T=n
  Y=rep(0,T)
  Y[1]=rbinom(1,size,prob)
  for(i in 1:(n-1))
  {
    u1=pbinom(Y[i],size,prob)
    u2=pbinom(Y[i]-1,size,prob)
    g=dbinom(Y[i],size,prob)
    
    repeat{
      u3=runif(1)
      s_function=function(y){ ((C0(y,u1,alpha)-C0(y,u2,alpha)))-g*u3 }
      ut=uniroot(s_function,lower=0,upper =1,extendInt = c("yes"))
      if(ut$root<1&&ut$root>0){break}
    }
    
    Y[i+1]=qbinom(ut$root,size,prob)
  }
  return(Y)

}