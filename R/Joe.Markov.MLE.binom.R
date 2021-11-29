Joe.Markov.MLE.binom = function (Y, size, k = 3, plot = TRUE, GOF = FALSE) 
{
  
  n = length(Y)
  
  log.l = function(par) {
    p = 1/(exp(par[1]) + 1)
    alpha = exp(par[2]) + 1
    
    A = function(u,v){
      (1-u)^alpha+(1-v)^alpha-(1-u)^alpha*(1-v)^alpha
    }
    
    l = log(dbinom(x = Y[1], size = size, prob = p)) - 
      sum(log(dbinom(x = Y[1:(n - 1)], size = size, prob = p))) + 
      sum(log(A(pbinom(q = Y[2:n], size = size, prob = p),
                pbinom(q = Y[1:(n-1)]-1, size = size, prob = p))^(1/alpha) - 
                A(pbinom(q = Y[2:n], size = size, prob = p),
                  pbinom(q = Y[1:(n-1)], size = size, prob = p))^(1/alpha) -
                A(pbinom(q = Y[2:n]-1, size = size, prob = p),
                  pbinom(q = Y[1:(n-1)]-1, size = size, prob = p))^(1/alpha) + 
                A(pbinom(q = Y[2:n]-1, size = size, prob = p),
                  pbinom(q = Y[1:(n-1)], size = size, prob = p))^(1/alpha)))
    return(-l)
  }
  alpha_0 = 2
  p_0 = mean(Y)/size
  res = nlm(f = log.l, p = c(log(1/p_0 - 1), log(alpha_0 + 1)), hessian = TRUE)
  p = 1/(exp(res$estimate[1]) + 1)
  alpha = exp(res$estimate[2]) + 1
  inverse_Hessian = solve(-res$hessian, tol = 10^(-50))
  SEP = sqrt(-inverse_Hessian[1, 1])
  SEp = SEP * ((p * (1 - p)))
  SEA = sqrt(-inverse_Hessian[2, 2])
  SEalpha = SEA * (alpha - 1)
  lower_p = 1/((((1/(p)) - 1) * exp(SEP * qnorm(0.975))) + 1)
  upper_p = 1/((((1/(p)) - 1) * exp(-SEP * qnorm(0.975))) + 1)
  lower_alpha = (alpha - 1) * exp(-SEA * qnorm(0.975)) + 1
  upper_alpha = (alpha - 1) * exp(SEA * qnorm(0.975)) + 1
  mu = size * p
  sigma = sqrt(size * p * (1 - p))
  result_p = c(estimate = p, SE = SEp, Lower = lower_p, Upper = upper_p)
  result_a = c(estimate = alpha, SE = SEalpha, Lower = lower_alpha, Upper = upper_alpha)
  result = c(mu = mu, sigma = sigma)
  
  ###
  names(mu) = NULL
  names(sigma) = NULL
  UCL = mu + k * sigma
  LCL = mu - k * sigma
  if (plot == TRUE) {
    Min = min(min(Y), LCL)
    Max = max(max(Y), UCL)
    ts.plot(Y, type = "b", ylab = "Y", ylim = c(Min, 
                                                Max))
    abline(h = mu)
    abline(h = UCL, lty = "dotted", lwd = 2)
    abline(h = LCL, lty = "dotted", lwd = 2)
    text(0, LCL + (mu - LCL) * 0.1, "LCL")
    text(0, UCL - (UCL - mu) * 0.1, "UCL")
  }
  out_control = which((Y < LCL) | (UCL < Y))
  if (length(out_control) == 0) {
    out_control = "NONE"
  }
  n = length(Y)
  F_par = pbinom(sort(Y), size = size, prob = result_p[1])
  F_emp = 1:n/n
  if (GOF == TRUE) {
    plot(F_emp, F_par, xlab = "F_empirical", ylab = "F_parametric", 
         xlim = c(0, 1), ylim = c(0, 1))
    lines(x = c(0, 1), y = c(0, 1))
  }
  CM.test = sum((F_emp - F_par)^2)
  KS.test = max(abs(F_emp - F_par))
  list(p = result_p, alpha = result_a, 
       Control_Limit = c(Center = mu, Lower = LCL, Upper = UCL), 
       out_of_control = out_control, 
       Gradient = res$gradient, Hessian = res$hessian, Eigenvalue_Hessian = eigen(res$hessian)$value, 
       KS.test = KS.test, CM.test = CM.test, log_likelihood = -res$minimum)
}

