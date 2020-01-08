Clayton.Markov2.MLE = function(Y, k = 3, D = 1, plot = TRUE, GOF=FALSE){

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
  initial = c(mean(Y),log(sd(Y)), log(ifelse(tau_0 < 0, 1, 2*tau_0/(1-tau_0))))

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


  CL = c(Center = mu.hat, Lower = LCL, Upper = UCL)
  
  Gradient = res$gradient
  Hessian = -res$hessian

  inverse_Hessian=solve(Hessian,tol=10^(-50))
  SE_mu = sqrt(-inverse_Hessian[1,1])
  SE_sigma = sqrt(-inverse_Hessian[2,2])*sigma.hat
  SE_alpha = sqrt(-inverse_Hessian[3,3])*alpha.hat

  lower_mu = mu.hat-1.96*SE_mu
  upper_mu = mu.hat+1.96*SE_mu

  lower_sigma=sigma.hat*exp(-1.96*SE_sigma/sigma.hat)
  upper_sigma=sigma.hat*exp(1.96*SE_sigma/sigma.hat)

  lower_alpha=alpha.hat*exp(-1.96*SE_alpha/alpha.hat)
  upper_alpha=alpha.hat*exp(1.96*SE_alpha/alpha.hat)

  result.mu = c(estimate = mu.hat, SE = SE_mu, Lower = lower_mu, Upper = upper_mu)
  result.sigma = c(estimate = sigma.hat, SE = SE_sigma, Lower = lower_sigma, Upper = upper_sigma)
  result.alpha = c(estimate = alpha.hat, SE = SE_alpha, Lower = lower_alpha, Upper = upper_alpha)

  ###plot
  if(plot==TRUE){
    par(mar=c(4,5,2,5))
    plot(Y~c(1:n), type = "b", ylim = c(1.1*LCL-0.1*UCL, 1.1*UCL-0.1*LCL),
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

  ### Goodness-of-fit ###
  F_par=pnorm( (sort(Y)-mu.hat)/sigma.hat )
  F_emp=1:n/n

  CM.test=sum( (F_emp-F_par)^2 )
  KS.test=max( abs( F_emp-F_par ) )

  if(GOF==TRUE){
    plot(F_emp,F_par,xlab="F_empirical",ylab="F_parametric",xlim=c(0,1),ylim=c(0,1))
    lines(x = c(0,1), y = c(0,1))
  }

  ##output
  return(list( mu=result.mu, sigma = result.sigma, alpha = result.alpha,
               Control_Limit = CL, out_of_control = OC,
               Gradient = Gradient, Hessian = Hessian, Eigenvalue_Hessian=eigen(Hessian)$value,
               CM.test=CM.test, KS.test=KS.test,log_likelihood = -logL(res$estimate)))
}
