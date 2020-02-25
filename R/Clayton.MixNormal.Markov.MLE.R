Clayton.MixNormal.Markov.MLE=function(y)
{
  #----likelihood function----
  logL = function(alpha,mu1,mu2,sigma1,sigma2,p){
    n = length(y)
    res = sum( log( f_MixN(x = y, mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, p = p) ) ) +
      + (n-1) * log(1+alpha) -
      (1+alpha) * sum( log( F_MixN(x = y[1:(n-1)], mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, p = p) ) ) -
      (1+alpha) * sum( log( F_MixN(x = y[2:n], mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, p = p) ) ) -
      (1/alpha+2) * sum( log(
        F_MixN(x = y[1:(n-1)], mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, p = p)^(-alpha) +
          F_MixN(x = y[2:n], mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, p = p)^(-alpha) - 1
      ) )
    return(res)
  }
  #----function----
  f_MixN=function(x,mu1,mu2,sigma1,sigma2,p) #Mixture Normal p.d.f
  {
    p*dnorm(x,mu1,sigma1) + (1-p)*dnorm(x,mu2,sigma2)
  }

  F_MixN=function(x,mu1,mu2,sigma1,sigma2,p) #Mixture Normal c.d.f
  {
    p*pnorm(x,mu1,sigma1) + (1-p)*pnorm(x,mu2,sigma2)
  }

  F_inverse=function(u,mu1,mu2,sigma1,sigma2,p) #Mixture Normal inverse c.d.f
  {
    Equ=function(z){F_MixN(z,mu1,mu2,sigma1,sigma2,p)-u}
    y=uniroot(Equ,c(-100,100))
    y$root
  }
  #----function of Mixture Normal---------------------
  fn=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    p*dnorm(x,m1,sigma1)+(1-p)*dnorm(x,m2,sigma2)}
  Fn=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    p*pnorm(x,m1,sigma1)+(1-p)*pnorm(x,m2,sigma2)}

  fn1=function(x,m1,s1) {
    sigma1=exp(s1)
    dnorm(x,m1,sigma1)}
  fn2=function(x,m2,s2) {
    sigma2=exp(s2)
    dnorm(x,m2,sigma2)}
  Fn1=function(x,m1,s1) {
    sigma1=exp(s1)
    pnorm(x,m1,sigma1)}
  Fn2=function(x,m2,s2) {
    sigma2=exp(s2)
    pnorm(x,m2,sigma2)}

  #----partial function of Mixture Normal----------------
  f_m1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    p*(x-m1)*fn1(x,m1,s1)/sigma1^2}
  f_m2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (1-p)*(x-m2)*fn2(x,m2,s2)/sigma2^2}
  f_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((-p/sigma1)*fn1(x,m1,s1) + p*(x-m1)^2*fn1(x,m1,s1)/sigma1^3)*exp(s1)}
  f_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((-(1-p)/sigma2)*fn2(x,m2,s2) + (1-p)*(x-m2)^2*fn1(x,m2,s2)/sigma2^3)*exp(s2)}
  f_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (fn1(x,m1,s1)-fn2(x,m2,s2))*(-exp(-exp(q)+q))}

  F_m1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    -p*fn1(x,m1,s1)}
  F_m2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    -(1-p)*fn2(x,m2,s2)}
  F_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-p*(x-m1)/sigma1*fn1(x,m1,s1))*exp(s1)}
  F_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-(1-p)*(x-m2)/sigma2*fn2(x,m2,s2))*exp(s2)}
  F_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (Fn1(x,m1,s1)-Fn2(x,m2,s2))*(-exp(-exp(q)+q))}

  #----2nd partial function of Mixture Normal------------
  f_m1_m1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    p*(x-m1)^2/sigma1^4*fn1(x,m1,s1)-p/sigma1^2*fn1(x,m1,s1)}
  f_m1_m2=function(x,m1,m2,s1,s2,q) {0}
  f_m1_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-3*p*(x-m1)/sigma1^3*fn1(x,m1,s1) + p*(x-m1)^3*fn1(x,m1,s1)/sigma1^5)*exp(s1)}
  f_m1_s2=function(x,m1,m2,s1,s2,q) {0}
  f_m1_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((x-m1)*fn1(x,m1,s1)/sigma1^2)*(-exp(-exp(q)+q))}

  f_m2_m2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (1-p)*(x-m2)^2/sigma2^4*fn2(x,m2,s2)-(1-p)/sigma2^2*fn2(x,m2,s2)}
  f_m2_s1=function(x,m1,m2,s1,s2,q) {0}
  f_m2_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-3*(1-p)*(x-m2)/sigma2^3*fn2(x,m2,s2) + (1-p)*(x-m2)^3*fn2(x,m2,s2)/sigma2^5)*exp(s2)}
  f_m2_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-(x-m2)*fn2(x,m2,s2)/sigma2^2)*(-exp(-exp(q)+q))}

  f_s1_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((2*p/sigma1^2 - 5*p*(x-m1)^2/sigma1^4 + p*(x-m1)^4/sigma1^6)*fn1(x,m1,s1))*exp(2*s1) + f_s1(x,m1,m2,s1,s2,q)}

  f_s1_s2=function(x,m1,m2,s1,s2,q) {0}
  f_s1_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((-1/sigma1)*fn1(x,m1,s1) + (x-m1)^2*fn1(x,m1,s1)/sigma1^3) *exp(s1) *(-exp(-exp(q)+q))}
  f_s2_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((2*(1-p)/sigma2^2 - 5*(1-p)*(x-m2)^2/sigma2^4 + (1-p)*(x-m2)^4/sigma2^6)*fn2(x,m2,s2))*exp(2*s2) + f_s2(x,m1,m2,s1,s2,q)}
  f_s2_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((1/sigma2)*fn2(x,m2,s2) - (x-m2)^2*fn2(x,m2,s2)/sigma2^3) *exp(s2) *(-exp(-exp(q)+q))}
  f_q_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    f_q(x,m1,m2,s1,s2,q)*(-exp(q)+1)}

  F_m1_m1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    -p*(x-m1)/sigma1^2*fn1(x,m1,s1)}
  F_m1_m2=function(x,m1,m2,s1,s2,q) {0}
  F_m1_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (p/sigma1*fn1(x,m1,s1)-p*(x-m1)^2/sigma1^3*fn1(x,m1,s1))*exp(s1)}
  F_m1_s2=function(x,m1,m2,s1,s2,q) {0}
  F_m1_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-fn1(x,m1,s1))*(-exp(-exp(q)+q))}

  #---
  F_m2_m2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    -(1-p)*(x-m2)/sigma2^2*fn2(x,m2,s2)}
  F_m2_s1=function(x,m1,m2,s1,s2,q) {0}
  F_m2_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((1-p)/sigma2*fn2(x,m2,s2)-(1-p)*(x-m2)^2/sigma2^3*fn2(x,m2,s2))*exp(s2)}
  F_m2_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    fn2(x,m2,s2)*(-exp(-exp(q)+q))}

  #---
  F_s1_s1=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((2*p*(x-m1)/sigma1^2 - p*(x-m1)^3/sigma1^4)*fn1(x,m1,s1))*exp(2*s1) + F_s1(x,m1,m2,s1,s2,q)}
  F_s1_s2=function(x,m1,m2,s1,s2,q) {0}
  F_s1_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    (-(x-m1)/sigma1*fn1(x,m1,s1))*exp(s1)*(-exp(-exp(q)+q))}
  #---
  F_s2_s2=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((2*(1-p)*(x-m2)/sigma2^2 - (1-p)*(x-m2)^3/sigma2^4)*fn2(x,m2,s2))*exp(2*s2) + F_s2(x,m1,m2,s1,s2,q)}

  F_s2_q=function(x,m1,m2,s1,s2,q) {
    p=exp(-exp(q))
    sigma1=exp(s1)
    sigma2=exp(s2)
    ((x-m2)/sigma2*fn2(x,m2,s2))*exp(s2)*(-exp(-exp(q)+q))}
  #---
  F_q_q=function(x,m1,m2,s1,s2,q) {F_q(x,m1,m2,s1,s2,q)*(-exp(q)+1)}

  #----partial function of Log-Likelihood---------
  vv1=function(y,a,m1,m2,s1,s2,q)  #a
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    A=(1/(1+alpha)-log(U1)-log(U2)+(1/alpha+2)*(U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2))/(U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    B=(log(U1^(-alpha)+U2^(-alpha)-1)/alpha^2)*exp(a)
    sum(A)+sum(B)
  }

  vv2=function(y,a,m1,m2,s1,s2,q)  #m1
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]

    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)

    A=f_m1(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)
    B=-(1+alpha)*(UU1/U1)-(1+alpha)*(UU2/U2)+(1+2*alpha)/(U1^(-alpha)+U2^(-alpha)-1)*(U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)
    sum(A)+sum(B)
  }
  vv3=function(y,a,m1,m2,s1,s2,q)  #m2
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    A=f_m2(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)

    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)

    B=-(1+alpha)*(UU1/U1)-(1+alpha)*(UU2/U2)+(1/alpha+2)*(alpha)/(U1^(-alpha)+U2^(-alpha)-1)*(U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)
    sum(A)+sum(B)
  }
  vv4=function(y,a,m1,m2,s1,s2,q)  #s1
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    A=f_s1(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)

    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s1(y1,m1,m2,s1,s2,q)
    UU2=F_s1(y2,m1,m2,s1,s2,q)

    B=-(1+alpha)*(UU1/U1)-(1+alpha)*(UU2/U2)+(1/alpha+2)*(alpha)/(U1^(-alpha)+U2^(-alpha)-1)*(U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)
    sum(A)+sum(B)
  }
  vv5=function(y,a,m1,m2,s1,s2,q)  #s2
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    A=f_s2(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)

    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s2(y1,m1,m2,s1,s2,q)
    UU2=F_s2(y2,m1,m2,s1,s2,q)

    B=-(1+alpha)*(UU1/U1)-(1+alpha)*(UU2/U2)+(1/alpha+2)*(alpha)/(U1^(-alpha)+U2^(-alpha)-1)*(U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)
    sum(A)+sum(B)
  }
  vv6=function(y,a,m1,m2,s1,s2,q)  #q
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    A=f_q(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)

    U1 =Fn (y1,m1,m2,s1,s2,q)
    U2 =Fn (y2,m1,m2,s1,s2,q)
    UU1=F_q(y1,m1,m2,s1,s2,q)
    UU2=F_q(y2,m1,m2,s1,s2,q)

    B=-(1+alpha)*(UU1/U1)-(1+alpha)*(UU2/U2)+(1/alpha+2)*(alpha)/(U1^(-alpha)+U2^(-alpha)-1)*(U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)
    sum(A)+sum(B)
  }


  #----2rd order partial function of Log-Likelihood
  hh11=function(y,a,m1,m2,s1,s2,q) #alpha_alpha
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    A= (-1/(1+alpha)^2 - 2/alpha^3*log(U1^(-alpha)+U2^(-alpha)-1) - (1/alpha^2)*(U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2))/(U1^(-alpha)+U2^(-alpha)-1))*exp(2*a)
    B= (-(1/alpha+2)*((U1^(-alpha)*log(U1)^2+U2^(-alpha)*log(U2)^2)*(U1^(-alpha)+U2^(-alpha)-1)-(U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2))^2)/(U1^(-alpha)+U2^(-alpha)-1)^2)*exp(2*a)
    C= (-(1/alpha^2)*(U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2))/(U1^(-alpha)+U2^(-alpha)-1))*exp(2*a)
    D= vv1(y,a,m1,m2,s1,s2,q)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh12=function(y,a,m1,m2,s1,s2,q) #alpha_m1  OK
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)

    A= (-UU1/U1-UU2/U2)*exp(a)
    B= (           2 * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) / (U1^(-alpha)+U2^(-alpha)-1)) *exp(a)
    C= (+(1+2*alpha) * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2)) / (U1^(-alpha)+U2^(-alpha)-1)^2)*exp(a)
    D= (-(1+2*alpha) * (U1^(-alpha-1)*UU1*log(U1)+U2^(-alpha-1)*UU2*log(U2))                             / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh13=function(y,a,m1,m2,s1,s2,q) #alpha_m2
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)

    A= (-UU1/U1-UU2/U2)*exp(a)
    B= (           2 * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    C= (+(1+2*alpha) * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2)) / (U1^(-alpha)+U2^(-alpha)-1)^2)*exp(a)
    D= (-(1+2*alpha) * (U1^(-alpha-1)*UU1*log(U1)+U2^(-alpha-1)*UU2*log(U2))                             / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh14=function(y,a,m1,m2,s1,s2,q) #alpha_s1
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s1(y1,m1,m2,s1,s2,q)
    UU2=F_s1(y2,m1,m2,s1,s2,q)

    A= (-UU1/U1-UU2/U2)*exp(a)
    B= (           2 * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    C= (+(1+2*alpha) * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2)) / (U1^(-alpha)+U2^(-alpha)-1)^2)*exp(a)
    D= (-(1+2*alpha) * (U1^(-alpha-1)*UU1*log(U1)+U2^(-alpha-1)*UU2*log(U2))                             / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh15=function(y,a,m1,m2,s1,s2,q) #alpha_s2
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s2(y1,m1,m2,s1,s2,q)
    UU2=F_s2(y2,m1,m2,s1,s2,q)

    A= (-UU1/U1-UU2/U2)*exp(a)
    B= (           2 * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    C= (+(1+2*alpha) * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2)) / (U1^(-alpha)+U2^(-alpha)-1)^2)*exp(a)
    D= (-(1+2*alpha) * (U1^(-alpha-1)*UU1*log(U1)+U2^(-alpha-1)*UU2*log(U2))                             / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh16=function(y,a,m1,m2,s1,s2,q) #alpha_q
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_q(y1,m1,m2,s1,s2,q)
    UU2=F_q(y2,m1,m2,s1,s2,q)

    A= (-UU1/U1-UU2/U2)*exp(a)
    B= (           2 * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    C= (+(1+2*alpha) * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha)*log(U1)+U2^(-alpha)*log(U2)) / (U1^(-alpha)+U2^(-alpha)-1)^2)*exp(a)
    D= (-(1+2*alpha) * (U1^(-alpha-1)*UU1*log(U1)+U2^(-alpha-1)*UU2*log(U2))                             / (U1^(-alpha)+U2^(-alpha)-1))*exp(a)
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh22=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)
    UUU1=F_m1_m1(y1,m1,m2,s1,s2,q)
    UUU2=F_m1_m1(y2,m1,m2,s1,s2,q)

    A= f_m1_m1(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)-f_m1(y,m1,m2,s1,s2,q)^2/fn(y,m1,m2,s1,s2,q)^2
    B= - (1+alpha)*(UUU1/U1-UU1^2/U1^2) - (1+alpha)*(UUU2/U2-UU2^2/U2^2)
    C= (1+2*alpha) * ((-alpha-1)*U1^(-alpha-2)*UU1^2+(-alpha-1)*U2^(-alpha-2)*UU2^2 + U1^(-alpha-1)*UUU1+U2^(-alpha-1)*UUU2) / (U1^(-alpha)+U2^(-alpha)-1)
    D= (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)^2 / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh23=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)
    uu1=F_m2(y1,m1,m2,s1,s2,q)
    uu2=F_m2(y2,m1,m2,s1,s2,q)

    A= -f_m1(y,m1,m2,s1,s2,q) * f_m2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B=  (1+alpha)*UU1*uu1/U1^2 + (1+alpha)*UU2*uu2/U2^2
    C= -(1+2*alpha) * (alpha+1) * (U1^(-alpha-2)*UU1*uu1+U2^(-alpha-2)*UU2*uu2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh24=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)
    uu1=F_s1(y1,m1,m2,s1,s2,q)
    uu2=F_s1(y2,m1,m2,s1,s2,q)
    XX1=F_m1_s1(y1,m1,m2,s1,s2,q)
    XX2=F_m1_s1(y2,m1,m2,s1,s2,q)

    A= f_m1_s1(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_m1(y,m1,m2,s1,s2,q) * f_s1(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh25=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)
    uu1=F_s2(y1,m1,m2,s1,s2,q)
    uu2=F_s2(y2,m1,m2,s1,s2,q)

    A= - f_m1(y,m1,m2,s1,s2,q) * f_s2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= + (1+alpha)*UU1*uu1/U1^2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh26=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m1(y1,m1,m2,s1,s2,q)
    UU2=F_m1(y2,m1,m2,s1,s2,q)
    uu1=F_q(y1,m1,m2,s1,s2,q)
    uu2=F_q(y2,m1,m2,s1,s2,q)
    XX1=F_m1_q(y1,m1,m2,s1,s2,q)
    XX2=F_m1_q(y2,m1,m2,s1,s2,q)

    A= f_m1_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_m1(y,m1,m2,s1,s2,q) * f_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh33=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)
    UUU1=F_m2_m2(y1,m1,m2,s1,s2,q)
    UUU2=F_m2_m2(y2,m1,m2,s1,s2,q)

    A= f_m2_m2(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)-f_m2(y,m1,m2,s1,s2,q)^2/fn(y,m1,m2,s1,s2,q)^2
    B= - (1+alpha)*(UUU1/U1-UU1^2/U1^2) - (1+alpha)*(UUU2/U2-UU2^2/U2^2)
    C= (1+2*alpha) * ((-alpha-1)*U1^(-alpha-2)*UU1^2+(-alpha-1)*U2^(-alpha-2)*UU2^2 + U1^(-alpha-1)*UUU1+U2^(-alpha-1)*UUU2) / (U1^(-alpha)+U2^(-alpha)-1)
    D= (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)^2 / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh34=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)
    uu1=F_s1(y1,m1,m2,s1,s2,q)
    uu2=F_s1(y2,m1,m2,s1,s2,q)
    XX1=F_m2_s1(y1,m1,m2,s1,s2,q)
    XX2=F_m2_s1(y2,m1,m2,s1,s2,q)

    A= f_m2_s1(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_m2(y,m1,m2,s1,s2,q) * f_s1(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh35=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)
    uu1=F_s2(y1,m1,m2,s1,s2,q)
    uu2=F_s2(y2,m1,m2,s1,s2,q)
    XX1=F_m2_s2(y1,m1,m2,s1,s2,q)
    XX2=F_m2_s2(y2,m1,m2,s1,s2,q)

    A= f_m2_s2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_m2(y,m1,m2,s1,s2,q) * f_s2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh36=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_m2(y1,m1,m2,s1,s2,q)
    UU2=F_m2(y2,m1,m2,s1,s2,q)
    uu1=F_q(y1,m1,m2,s1,s2,q)
    uu2=F_q(y2,m1,m2,s1,s2,q)
    XX1=F_m2_q(y1,m1,m2,s1,s2,q)
    XX2=F_m2_q(y2,m1,m2,s1,s2,q)

    A= f_m2_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_m2(y,m1,m2,s1,s2,q) * f_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh44=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s1(y1,m1,m2,s1,s2,q)
    UU2=F_s1(y2,m1,m2,s1,s2,q)
    UUU1=F_s1_s1(y1,m1,m2,s1,s2,q)
    UUU2=F_s1_s1(y2,m1,m2,s1,s2,q)

    A= f_s1_s1(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)-f_s1(y,m1,m2,s1,s2,q)^2/fn(y,m1,m2,s1,s2,q)^2
    B= - (1+alpha)*(UUU1/U1-UU1^2/U1^2) - (1+alpha)*(UUU2/U2-UU2^2/U2^2)
    C= (1+2*alpha) * ((-alpha-1)*U1^(-alpha-2)*UU1^2+(-alpha-1)*U2^(-alpha-2)*UU2^2 + U1^(-alpha-1)*UUU1+U2^(-alpha-1)*UUU2) / (U1^(-alpha)+U2^(-alpha)-1)
    D= (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)^2 / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh45=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s1(y1,m1,m2,s1,s2,q)
    UU2=F_s1(y2,m1,m2,s1,s2,q)
    uu1=F_s2(y1,m1,m2,s1,s2,q)
    uu2=F_s2(y2,m1,m2,s1,s2,q)
    XX1=F_s1_s2(y1,m1,m2,s1,s2,q)
    XX2=F_s1_s2(y2,m1,m2,s1,s2,q)

    A= f_s1_s2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_s1(y,m1,m2,s1,s2,q) * f_s2(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh46=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s1(y1,m1,m2,s1,s2,q)
    UU2=F_s1(y2,m1,m2,s1,s2,q)
    uu1=F_q(y1,m1,m2,s1,s2,q)
    uu2=F_q(y2,m1,m2,s1,s2,q)
    XX1=F_s1_q(y1,m1,m2,s1,s2,q)
    XX2=F_s1_q(y2,m1,m2,s1,s2,q)

    A= f_s1_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_s1(y,m1,m2,s1,s2,q) * f_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh55=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s2(y1,m1,m2,s1,s2,q)
    UU2=F_s2(y2,m1,m2,s1,s2,q)
    UUU1=F_s2_s2(y1,m1,m2,s1,s2,q)
    UUU2=F_s2_s2(y2,m1,m2,s1,s2,q)

    A= f_s2_s2(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)-f_s2(y,m1,m2,s1,s2,q)^2/fn(y,m1,m2,s1,s2,q)^2
    B= - (1+alpha)*(UUU1/U1-UU1^2/U1^2) - (1+alpha)*(UUU2/U2-UU2^2/U2^2)
    C= (1+2*alpha) * ((-alpha-1)*U1^(-alpha-2)*UU1^2+(-alpha-1)*U2^(-alpha-2)*UU2^2 + U1^(-alpha-1)*UUU1+U2^(-alpha-1)*UUU2) / (U1^(-alpha)+U2^(-alpha)-1)
    D= (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)^2 / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh56=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_s2(y1,m1,m2,s1,s2,q)
    UU2=F_s2(y2,m1,m2,s1,s2,q)
    uu1=F_q(y1,m1,m2,s1,s2,q)
    uu2=F_q(y2,m1,m2,s1,s2,q)
    XX1=F_s2_q(y1,m1,m2,s1,s2,q)
    XX2=F_s2_q(y2,m1,m2,s1,s2,q)

    A= f_s2_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q) - f_s2(y,m1,m2,s1,s2,q) * f_q(y,m1,m2,s1,s2,q) / fn(y,m1,m2,s1,s2,q)^2
    B= -(1+alpha)*XX1/U1 + (1+alpha)*UU1*uu1/U1^2 -(1+alpha)*XX2/U2 + (1+alpha)*UU2*uu2/U2^2
    C=  (1+2*alpha) *  ((-alpha-1)*U1^(-alpha-2)*UU1*uu1 + (-alpha-1)*U2^(-alpha-2)*UU2*uu2 + U1^(-alpha-1)*XX1 + U2^(-alpha-1)*XX2) / (U1^(-alpha)+U2^(-alpha)-1)
    D=  (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2) * (U1^(-alpha-1)*uu1+U2^(-alpha-1)*uu2) / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  hh66=function(y,a,m1,m2,s1,s2,q)
  {
    alpha=exp(a)-1
    y1=y[-length(y)]
    y2=y[-1]
    U1=Fn(y1,m1,m2,s1,s2,q)
    U2=Fn(y2,m1,m2,s1,s2,q)
    UU1=F_q(y1,m1,m2,s1,s2,q)
    UU2=F_q(y2,m1,m2,s1,s2,q)
    UUU1=F_q_q(y1,m1,m2,s1,s2,q)
    UUU2=F_q_q(y2,m1,m2,s1,s2,q)

    A= f_q_q(y,m1,m2,s1,s2,q)/fn(y,m1,m2,s1,s2,q)-f_q(y,m1,m2,s1,s2,q)^2/fn(y,m1,m2,s1,s2,q)^2
    B= - (1+alpha)*(UUU1/U1-UU1^2/U1^2) - (1+alpha)*(UUU2/U2-UU2^2/U2^2)
    C= (1+2*alpha) * ((-alpha-1)*U1^(-alpha-2)*UU1^2+(-alpha-1)*U2^(-alpha-2)*UU2^2 + U1^(-alpha-1)*UUU1+U2^(-alpha-1)*UUU2) / (U1^(-alpha)+U2^(-alpha)-1)
    D= (1+2*alpha) * alpha * (U1^(-alpha-1)*UU1+U2^(-alpha-1)*UU2)^2 / (U1^(-alpha)+U2^(-alpha)-1)^2
    sum(A)+sum(B)+sum(C)+sum(D)
  }

  #----vector and Hessian matrix----
  vv=function(y,a,m1,m2,s1,s2,q)
  {
    vec=matrix(0,6,1)
    vec[1,1]=vv1(y,a,m1,m2,s1,s2,q)
    vec[2,1]=vv2(y,a,m1,m2,s1,s2,q)
    vec[3,1]=vv3(y,a,m1,m2,s1,s2,q)
    vec[4,1]=vv4(y,a,m1,m2,s1,s2,q)
    vec[5,1]=vv5(y,a,m1,m2,s1,s2,q)
    vec[6,1]=vv6(y,a,m1,m2,s1,s2,q)

    return(vec)
  }

  hh=function(y,a,m1,m2,s1,s2,q)
  {
    hes=matrix(0,6,6)
    hes[1,1]=hh11(y,a,m1,m2,s1,s2,q)

    hes[2,1]=hh12(y,a,m1,m2,s1,s2,q)
    hes[1,2]=hh12(y,a,m1,m2,s1,s2,q)

    hes[2,2]=hh22(y,a,m1,m2,s1,s2,q)

    hes[3,1]=hh13(y,a,m1,m2,s1,s2,q)
    hes[1,3]=hh13(y,a,m1,m2,s1,s2,q)

    hes[2,3]=hh23(y,a,m1,m2,s1,s2,q)
    hes[3,2]=hh23(y,a,m1,m2,s1,s2,q)

    hes[4,1]=hh14(y,a,m1,m2,s1,s2,q)
    hes[1,4]=hh14(y,a,m1,m2,s1,s2,q)

    hes[5,1]=hh15(y,a,m1,m2,s1,s2,q)
    hes[1,5]=hh15(y,a,m1,m2,s1,s2,q)

    hes[2,4]=hh24(y,a,m1,m2,s1,s2,q)
    hes[4,2]=hh24(y,a,m1,m2,s1,s2,q)

    hes[3,3]=hh33(y,a,m1,m2,s1,s2,q)

    hes[6,1]=hh16(y,a,m1,m2,s1,s2,q)
    hes[1,6]=hh16(y,a,m1,m2,s1,s2,q)

    hes[5,2]=hh25(y,a,m1,m2,s1,s2,q)
    hes[2,5]=hh25(y,a,m1,m2,s1,s2,q)

    hes[3,4]=hh34(y,a,m1,m2,s1,s2,q)
    hes[4,3]=hh34(y,a,m1,m2,s1,s2,q)

    hes[2,6]=hh26(y,a,m1,m2,s1,s2,q)
    hes[6,2]=hh26(y,a,m1,m2,s1,s2,q)

    hes[3,5]=hh35(y,a,m1,m2,s1,s2,q)
    hes[5,3]=hh35(y,a,m1,m2,s1,s2,q)

    hes[4,4]=hh44(y,a,m1,m2,s1,s2,q)

    hes[3,6]=hh36(y,a,m1,m2,s1,s2,q)
    hes[6,3]=hh36(y,a,m1,m2,s1,s2,q)

    hes[4,5]=hh45(y,a,m1,m2,s1,s2,q)
    hes[5,4]=hh45(y,a,m1,m2,s1,s2,q)

    hes[4,6]=hh46(y,a,m1,m2,s1,s2,q)
    hes[6,4]=hh46(y,a,m1,m2,s1,s2,q)

    hes[5,5]=hh55(y,a,m1,m2,s1,s2,q)

    hes[5,6]=hh56(y,a,m1,m2,s1,s2,q)
    hes[6,5]=hh56(y,a,m1,m2,s1,s2,q)

    hes[6,6]=hh66(y,a,m1,m2,s1,s2,q)

    return(hes)
  }
  #----Estimate Method----
  Newton_Raphson=function(w)
  {
    random=function(x){runif(6,-0.1,0.1)}
    theta=matrix(w,6,1) #theta=c(a,mu1,mu2,s1,s2,q)
    output=c()
    error=0
    qq=0
    repeat
    {
      qq=qq+1
      theta_l = theta
      Q=
      {
        U1=Fn(y[-1]        ,theta[2,1],theta[3,1],theta[4,1],theta[5,1],theta[6,1])
        U2=Fn(y[-length(y)],theta[2,1],theta[3,1],theta[4,1],theta[5,1],theta[6,1])
        alpha=exp(theta[1,1])-1
        (U1^(-alpha)+U2^(-alpha)-1)
      }
      if(length(Q[Q<0])>0)
      {
        qq=0
        cat("e1 ")
        X=w+runif(6,-0.05,0.05)
        theta=matrix(X,6,1)
        theta_l = theta
        next
      }

      vector  = vv(y,theta[1,1],theta[2,1],theta[3,1],theta[4,1],theta[5,1],theta[6,1])
      hessian = hh(y,theta[1,1],theta[2,1],theta[3,1],theta[4,1],theta[5,1],theta[6,1])

      if(is.na(abs(det(hessian))) || abs(det(hessian))<1e-100 || abs(det(hessian))>1e+100 ||qq>100) #|| theta[1,1]>log(20)abs(exp(-exp(theta[6]))-0.5)>0.45)
      {
        qq=0
        cat("e2 ")
        X=w+runif(6,-0.05,0.05)
        theta=matrix(X,6,1)
        theta_l = theta
        next
      }

      theta  = theta - Matrix_I(hessian) %*% vector

      if(is.na(sum(theta)) || sum(theta_l)==0)
      {
        X=w+c(runif(1,-0.1,0.1),runif(4,-0.5,0.5),runif(1,-1,1))
        theta=matrix(X,6,1)
      }

      if(abs(sum(theta-theta_l))<1e-8)
      {
        z=theta
        eig=eigen(hh(y,z[1],z[2],z[3],z[4],z[5],z[6]))$values
        if(sum(sign(eig))==-6)
        {
          output=theta
          break
        }
        X=w+runif(6,-0.05,0.05)
        theta=matrix(X,6,1)
      }
    }
    return(output)
  }
  #-----k-means clustering----
  k_means_clustering_EM=function(y)
  {
    y
    z1=0.5*(min(y)+max(y))+sd(y)
    z2=0.5*(min(y)+max(y))-sd(y)
    z3=sd(y)*0.5
    z4=sd(y)*0.5
    z5=0.5
    #---
    t_=c()
    t=c(z1,z2,z3,z4,z5)
    r=c()
    x=0
    repeat
    {
      x=x+1
      t_=t
      r=(t[5]*dnorm(y,t[1],t[3]))/(t[5]*dnorm(y,t[1],t[3])+(1-t[5])*dnorm(y,t[2],t[4]))

      t[1]=sum(r*y)/sum(r)
      t[2]=sum((1-r)*y)/sum(1-r)
      t[3]=sqrt(sum(r*(y-t[1])^2)/sum(r))
      t[4]=sqrt(sum((1-r)*(y-t[2])^2)/sum(1-r))
      t[5]=sum(r)/length(r)

      if(sum(abs(t-t_))<1e-12 && t[3]!=0 && t[4]!=0)
      {
        if(t[1]>t[2])
        {
          s=t
          t=c(s[2],s[1],s[4],s[3],1-s[5])
        }
        break
      }
    }
    return(t)
  }

  Matrix_I=function(H)
  {
    l=sqrt(length(H))
    HH=matrix(0,l,l)
    D=det(H)
    for(i in 1:l)
    {
      for(j in 1:l)
      {
        HH[i,j]=det(H[-i,-j])*(-1)^(i+j)
      }
    }
    HH=t(HH/det(H))
    return(HH)
  }

  transform=function(t)
  {
    c(exp(t[1])-1,t[2],t[3],exp(t[4]),exp(t[5]),exp(-exp(t[6])))
  }

  transform_1=function(t)
  {
    c(log(t[1]+1),t[2],t[3],log(t[4]),log(t[5]),log(-log(t[6])))
  }
  #-------Estimating-----
  #---Estimate alpha----------------------------------------------
  alpha_hat=
  {
    M=matrix(0,length(y)-1,2)
    M[,1]=y[-length(y)]
    M[,2]=y[-1]
    M=M[order(M[,1]),]
    e=M[,2]

    cc=0
    for(i in 1:(length(e)-1))
    {
      ee=e[(i+1):length(e)]
      cc=cc+sum(sign(ee-e[i]))
    }
    tau=cc/(length(e)*(length(e)-1)/2)
    2*tau/(1-tau)
  }
  set.seed(NULL)
  #---k-means clustering-----------------------------------
  tt=k_means_clustering_EM(y)
  #---Newton Raphson for mixnormal---------------------------------
  o=c()
  w=transform_1(c(alpha_hat,tt))
  o=Newton_Raphson(w)
  if(o[2]>o[3])
  {
    oo=o
    o=c(oo[1],oo[3],oo[2],oo[5],oo[4],log(-log(1-exp(-exp(o[6])))))
  }
  tr=transform(o)
  o=transform_1(tr)
  #-----Standard Error-----
  #-----coverage-----
  H=Matrix_I(-hh(y,o[1],o[2],o[3],o[4],o[5],o[6])) #Information matrix
  ssd=c()
  for(t in 1:6)
  {
    if(t==1)#alpha
    {
      ssd[t]=(tr[t]+1)*sqrt(H[t,t])
    }
    if(t==2 || t==3)#mu*
    {
      ssd[t]=sqrt(H[t,t])
    }
    if(t==4 || t==5 )#sigma*
    {
      ssd[t]=tr[t]*sqrt(H[t,t])
    }
    if(t==6)#p
    {
      ssd[t]=-(tr[t]*log(tr[t]))*sqrt(H[t,t])
    }

  }

  result=c(alpha=tr[1],mu1=tr[2],mu2=tr[3],sigma1=tr[4],sigma2=tr[5],p=tr[6])
  Stand_err=c(alpha=ssd[1],mu1=ssd[2],mu2=ssd[3],sigma1=ssd[4],sigma2=ssd[5],p=ssd[6])
  gradient=vv(y,o[1],o[2],o[3],o[4],o[5],o[6])
  hessian=hh(y,o[1],o[2],o[3],o[4],o[5],o[6])
  eigen_hessian=eigen(hessian)$value

  inv_hessian = solve(hessian, tol = 10^(-50))

  res_alpha = c(estimate = tr[1], SE = ssd[1], Lower = tr[1]-1.96*ssd[1], Upper = tr[1]+1.96*ssd[1])
  res_mu1 = c(estimate = tr[2], SE = ssd[2], Lower = tr[2]-1.96*ssd[2], Upper = tr[2]+1.96*ssd[2])
  res_mu2 = c(estimate = tr[3], SE = ssd[3], Lower = tr[3]-1.96*ssd[3], Upper = tr[3]+1.96*ssd[3])
  res_sigma1 = c(estimate = tr[4], SE = ssd[4], Lower = tr[4]-1.96*ssd[4], Upper = tr[4]+1.96*ssd[4])
  res_sigma2 = c(estimate = tr[5], SE = ssd[5], Lower = tr[5]-1.96*ssd[5], Upper = tr[5]+1.96*ssd[5])
  res_p = c(estimate = tr[6], SE = ssd[6], Lower = tr[6]-1.96*ssd[6], Upper = tr[6]+1.96*ssd[6])

  return(
    list(alpha = res_alpha, mu1 = res_mu1, mu2 = res_mu2,
         sigma1 = res_sigma1, sigma2 = res_sigma2, p = res_p,
         Gradient=as.vector(gradient),Hessian=hessian,
         Eigenvalue_Hessian=eigen_hessian,
         log_likelihood = as.numeric(logL(alpha = result[1],mu1 = result[2], mu2 = result[3],
                               sigma1 = result[4], sigma2 = result[5], p = result[6])))
  )
}
