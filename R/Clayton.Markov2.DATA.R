Clayton.Markov2.DATA = function(n,mu,sigma,alpha){
  U = numeric(n)
  U[1] = runif(1, min = 0, max = 1)
  U[2] = ((runif(1, min = 0, max = 1)^(-alpha/(1+alpha))-1)*U[1]^(-alpha)+1)^(-1/alpha)
  for(i in c(3:n)){
    U[i] = ( runif(1, min = 0, max = 1)^(-alpha/(1+2*alpha))*(U[i-1]^(-alpha)+U[i-2]^(-alpha)-1)-
               U[i-1]^(-alpha)-U[i-2]^(-alpha)+2 )^(-1/alpha)
  }
  Y = qnorm(p = U, mean = mu, sd = sigma)
  return(Y)
}