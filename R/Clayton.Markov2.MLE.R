Clayton.Markov2.MLE = function(Y, k = 3, D = 1, plot = TRUE){
  
  n = length(Y)
  rec = NA
  
  ###log-likelihood
  logL = function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    alpha = exp(par[3])
    
    Z = (Y - mu)/sigma
    G.yt = pnorm(q = Z, mean = 0, sd = 1)
    g.ty = 1/sigma * dnorm(x = Z, mean = 0, sd = 1)
    
    l = (n-2)*log(1+2*alpha) + log(1+alpha) -
      (1/alpha+3)*sum(log(G.yt[3:n]^(-alpha)+G.yt[2:(n-1)]^(-alpha)+G.yt[1:(n-2)]^(-alpha)-2)) +
      (1/alpha+2)*sum(log(G.yt[3:(n-1)]^(-alpha)+G.yt[2:(n-2)]^(-alpha)-1)) -
      (alpha+1)*sum(log(G.yt)) + sum(log(g.ty))
    
    return(-l)
    
  }
  
  ###initial value and randomize
  tau_0 = cor(Y[2:n],Y[1:(n-1)],method="kendall")
  initial = c(mean(Y),log(sd(Y)), log(2*tau_0/(1-tau_0)))
  
  count = 0
  repeat{
    count = count + 1
    res = try(nlm(logL, initial , hessian = TRUE))
    if( class(res)!="try-error" ){
      break;
    }else{
      initial = initial + runif(n = 3, min = -D, max = D)
    }
    if(count>100){
      return(warning("error"))
      break;
    }
  }
  
  ###result
  mu.hat = res$estimate[1]
  sigma.hat = exp(res$estimate[2])
  alpha.hat = exp(res$estimate[3])
  UCL = mu.hat+k*sigma.hat
  LCL = mu.hat-k*sigma.hat
  
  ###plot
  if(plot==TRUE){
    par(mar=c(4,5,2,5))
    plot(Y~c(1:n), type = "b", ylim = c(1.1*UCL-0.1*LCL, 1.1*LCL-0.1*UCL),
         ylab = "Y", xlab = "Time", cex = 1, cex.lab = 1)
    abline(h=UCL, lty = 3, lwd = 2)
    abline(h=LCL, lty = 3, lwd = 2)
    abline(h=mu.hat, lty = 1, lwd = 1)
    axis(4, at = c(UCL, LCL, mu.hat), 
         labels = c("UCL", "LCL", "mu"), cex = 1, las = 1)
  }
  
  ###out of control
  OC = which((Y < LCL) | (UCL < Y))
  if (length(OC) == 0) {
    OC = "NONE"
  }
  
  ##output
  MLE = c(mu.hat, sigma.hat, alpha.hat, UCL, LCL)
  names(MLE) = c("mu", "sigma", "alpha", "UCL", "LCL")
  
  return(list( estimate = MLE, out_of_control = OC, gradient = res$gradient, hessian = res$hessian))
}