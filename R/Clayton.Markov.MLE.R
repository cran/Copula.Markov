Clayton.Markov.MLE <-
  function(Y,k=3,D=1,plot=TRUE,GOF=FALSE, method = "nlm"){
    n=length(Y)   ##sample size##
    log.L = function(par){

      mu = par[1]
      sigma = exp(par[2])
      alpha = exp(par[3])

      g_yt = dnorm(Y, mean = mu, sd = sigma)
      G_yt = pnorm(Y, mean = mu, sd = sigma)

      C_11 = function(u,v){
        (1+alpha)*(u^(-alpha)+v^(-alpha)-1)^(-1/alpha-2)*v^(-alpha-1)*u^(-alpha-1)
      }

      res = sum( log(g_yt) ) + sum( log( C_11( G_yt[-n], G_yt[-1] ) ) )

      return(-res)
    }

    if (method == "Newton"){
      G=function(y,mu,sigma){pnorm((y-mu)/sigma)}   #G function
      g=function(y,mu,sigma){dnorm((y-mu)/sigma)}   #g function
      ##############Log-likelihood function #################
      L_function=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        A=U_t_1^(-alpha)+U_t^(-alpha)-1
        Z=log(1+alpha)-(1+alpha)*log(U_t_1)-(1+alpha)*log(U_t)-(1/alpha+2)*log(A)
        ZZ=log(g(Y[1:n],mu,sigma)/sigma)
        return((sum(Z)+sum(ZZ))/n)
      }
      ##############dL/dmu############################
      dL_dmu=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        A1=u_t_1/U_t_1+u_t/U_t
        A2=(u_t_1*U_t_1^(-(alpha+1))+u_t*U_t^(-(alpha+1)))/(U_t_1^(-alpha)+U_t^(-alpha)-1)
        A=(alpha+1)/sigma*A1-(2*alpha+1)/sigma*A2
        A3=(Y[1:n]-mu)/sigma^2
        return((sum(A)+sum(A3))/n)
      }
      ############dL/dsigma###########################
      dL_dsigma=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        B1=(Y[1:n-1]-mu)/sigma^2*u_t_1/U_t_1
        B2=(Y[2:n]-mu)/sigma^2*u_t/U_t
        B3=(Y[1:n-1]-mu)/sigma^2*U_t_1^(-(alpha+1))*u_t_1
        B4=(Y[2:n]-mu)/sigma^2*U_t^(-(alpha+1))*u_t
        B5=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        B=(alpha+1)*B1+(alpha+1)*B2-(2*alpha+1)*(B3+B4)/B5
        B6=(Y[1:n]-mu)^2/sigma^3-1/sigma
        return((sum(B)+sum(B6))/n)
      }
      ############dL/dalpha#############################
      dL_dalpha=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        E1=log(U_t_1*U_t)
        E2=log(  U_t_1^(-alpha)+U_t^(-alpha)-1  )/alpha^2
        E3=U_t_1^(-alpha)*log(U_t_1)+U_t^(-alpha)*log(U_t)
        E5=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        E=1/(1+alpha)-E1+E2+(2+1/alpha)*E3/E5
        return(sum(E)/n)
      }
      ############F function##############################
      F=function(mu,sigma,alpha){
        c(dL_dmu(mu,sigma,alpha),dL_dsigma(mu,sigma,alpha),dL_dalpha(mu,sigma,alpha))
      }
      ############d^2L/dmu^2#################################
      d2L_dmu2=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        H1=((Y[1:n-1]-mu)/sigma^2*u_t_1*U_t_1+u_t_1^2/sigma)/U_t_1^2
        H2=((Y[2:n]-mu)/sigma^2*u_t*U_t+u_t^2/sigma)/U_t^2
        H3=U_t_1^(-(2+alpha))*u_t_1^2
        H4=(Y[1:n-1]-mu)/sigma^2*u_t_1*U_t_1^(-(1+alpha))
        H5=U_t^(-(2+alpha))*u_t^2
        H6=(Y[2:n]-mu)/sigma^2*u_t*U_t^(-(1+alpha))
        H7=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        H8=U_t_1^(-(1+alpha))*u_t_1+U_t^(-(1+alpha))*u_t
        H=(alpha+1)/sigma*(H1+H2)-(2*alpha+1)/sigma*(((alpha+1)/sigma*H3+H4+(alpha+1)/sigma*H5+H6)*H7-alpha/sigma*H8^2)/H7^2
        return((sum(H)+(-n/sigma^2))/n)
      }
      ############d^2L/dsigma^2#################################
      J=J1=J2=J3=J4=J5=J6=J7=J8=J9=J10=J11=J12=J13=J14=J15=J16=c()
      d2L_dsigma2=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        J1=(Y[1:n-1]-mu)/sigma^3*u_t_1/U_t_1
        J2=(Y[1:n-1]-mu)^2/sigma^2+(Y[1:n-1]-mu)/sigma*u_t_1/U_t_1
        J4=(Y[2:n]-mu)/sigma^3*u_t/U_t
        J5=(Y[2:n]-mu)^2/sigma^2+(Y[2:n]-mu)/sigma*u_t/U_t
        J7=(Y[1:n-1]-mu)/sigma^3*U_t_1^(-(alpha+1))*u_t_1
        J8=(Y[1:n-1]-mu)/sigma*u_t_1/U_t_1
        J9=(Y[1:n-1]-mu)^2/sigma^2
        J10=(Y[2:n]-mu)/sigma^3*U_t^(-(alpha+1))*u_t
        J11=(Y[2:n]-mu)/sigma*u_t/U_t
        J12=(Y[2:n]-mu)^2/sigma^2
        J13=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        J15=(Y[1:n-1]-mu)/sigma^2*U_t_1^(-(alpha+1))*u_t_1
        J16=(Y[2:n]-mu)/sigma^2*U_t^(-(alpha+1))*u_t
        J=(alpha+1)*(J1*(-2+J2)+J4*(-2+J5))-(2*alpha+1)*(  (J7*(-2+(alpha+1)*J8+J9)+J10*(-2+(alpha+1)*J11+J12))*J13-alpha*(J15+J16)^2  )/J13^2
        J14=-3*(Y[1:n]-mu)^2/sigma^4+1/sigma^2
        return((sum(J)+sum(J14))/n)
      }
      ############d^2L/dalpha^2#################################
      d2L_dalpha2=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        K1=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        K2=U_t_1^(-alpha)*log(U_t_1)+U_t^(-alpha)*log(U_t)
        K4=U_t_1^(-alpha)*log(U_t_1)^2+U_t^(-alpha)*log(U_t)^2
        K=-1/(1+alpha)^2-2/alpha^3*log(K1)-2/alpha^2*K2/K1+(1/alpha+2)*( K2^2/K1^2-K4/K1 )
        return(sum(K)/n)
      }
      ############d^2L/dmudsigma#################################
      d2L_dmudsigma=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        L1=-(1+alpha)/sigma^2*(u_t_1/U_t_1+u_t/U_t)
        L2=((Y[1:n-1]-mu)^2/sigma^3*u_t_1*U_t_1+(Y[1:n-1]-mu)/sigma^2*(u_t_1)^2)/U_t_1^2
        L3=((Y[2:n]-mu)^2/sigma^3*u_t*U_t+(Y[2:n]-mu)/sigma^2*(u_t)^2)/U_t^2
        L4=(U_t_1^(-1-alpha)*u_t_1+U_t^(-1-alpha)*u_t)/(U_t_1^(-alpha)+U_t^(-alpha)-1)
        L5=1/(U_t_1^(-alpha)+U_t^(-alpha)-1)^2
        L6=(alpha+1)*U_t_1^(-alpha-2)*(Y[1:n-1]-mu)/sigma^2*u_t_1^2
        L7=U_t_1^(-alpha-1)*u_t_1*(Y[1:n-1]-mu)^2/sigma^3
        L8=(alpha+1)*U_t^(-alpha-2)*(Y[2:n]-mu)/sigma^2*u_t^2
        L9=U_t^(-alpha-1)*u_t*(Y[2:n]-mu)^2/sigma^3
        L10=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        L11=U_t_1^(-1-alpha)*u_t_1+U_t^(-1-alpha)*u_t
        L12=(Y[1:n-1]-mu)/sigma^2*U_t_1^(-1-alpha)*u_t_1+(Y[2:n]-mu)/sigma^2*U_t^(-1-alpha)*u_t
        L=L1+(alpha+1)/sigma*(L2+L3)+(2*alpha+1)/sigma^2*L4-(2*alpha+1)/sigma*L5*( (L6+L7+L8+L9)*L10-alpha*L11*L12 )
        LL=-2*(Y[1:n]-mu)/sigma^3
        return((sum(L)+sum(LL))/n)
      }
      ############d^2L/dmudalpha#################################
      d2L_dmudalpha=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        M1=(u_t_1/U_t_1+u_t/U_t)/sigma
        M2=(U_t_1^(-(alpha+1))*u_t_1+U_t^(-(alpha+1))*u_t)/(U_t_1^(-alpha)+U_t^(-alpha)-1)*2/sigma
        M3=-U_t_1^(-(1+alpha))*log(U_t_1)*u_t_1-U_t^(-(1+alpha))*log(U_t)*u_t
        M4=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        M5=U_t_1^(-(alpha+1))*u_t_1+U_t^(-(alpha+1))*u_t
        M6=-U_t_1^(-alpha)*log(U_t_1)-U_t^(-alpha)*log(U_t)
        M=M1-M2-(2*alpha+1)/sigma*( M3*M4-M5*M6 )/M4^2
        return(sum(M)/n)
      }
      ############d^2L/dsigmadalpha#################################
      d2L_dsigmadalpha=function(mu,sigma,alpha){
        U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
        u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
        O1=(Y[1:n-1]-mu)/sigma^2*u_t_1/U_t_1+(Y[2:n]-mu)/sigma^2*u_t/U_t
        O3=(Y[1:n-1]-mu)/sigma^2*U_t_1^(-(1+alpha))*u_t_1
        O4=(Y[2:n]-mu)/sigma^2*U_t^(-(1+alpha))*u_t
        O5=(U_t_1^(-alpha)+U_t^(-alpha)-1)
        O6=U_t_1^(-alpha)*log(U_t_1)+U_t^(-alpha)*log(U_t)
        O=O1-2*(O3+O4)/O5-(2*alpha+1)*( (-log(U_t_1)*O3-log(U_t)*O4)*O5+(O3+O4)*O6 )/O5^2
        return(sum(O)/n)
      }
      #############Jacobian function######################################
      Ja=function(mu,sigma,alpha){
        AA=matrix( c(d2L_dmu2(mu,sigma,alpha),d2L_dmudsigma(mu,sigma,alpha),d2L_dmudalpha(mu,sigma,alpha),d2L_dmudsigma(mu,sigma,alpha),d2L_dsigma2(mu,sigma,alpha),d2L_dsigmadalpha(mu,sigma,alpha),d2L_dmudalpha(mu,sigma,alpha),d2L_dsigmadalpha(mu,sigma,alpha),d2L_dalpha2(mu,sigma,alpha)),3,3 )
        return(AA)
      }
      ###############Multivariate Newton Raphson#####################
      X=matrix(,1,3)
      tau=cor(Y[1:n-1],Y[2:n],method="kendall")
      alpha_est=-2*tau/(tau-1)
      X[1,]=c(mean(Y),sd(Y),alpha_est)   #initial value
      i=2
      Ran.num=1
      iter.num = 0
      repeat{
        Z=X
        X=matrix(,i,3)
        X[1:i-1,]=Z[1:i-1,]
        ##
        Aa=Ja(X[i-1,1],X[i-1,2],X[i-1,3])
        Ainv11=Aa[2,2]*Aa[3,3]-Aa[3,2]*Aa[2,3]
        Ainv12=Aa[1,2]*Aa[3,3]-Aa[3,2]*Aa[1,3]
        Ainv13=Aa[1,2]*Aa[2,3]-Aa[2,2]*Aa[1,3]
        Ainv21=Aa[2,1]*Aa[3,3]-Aa[3,1]*Aa[2,3]
        Ainv22=Aa[1,1]*Aa[3,3]-Aa[3,1]*Aa[1,3]
        Ainv23=Aa[1,1]*Aa[2,3]-Aa[1,3]*Aa[2,1]
        Ainv31=Aa[2,1]*Aa[3,2]-Aa[3,1]*Aa[2,2]
        Ainv32=Aa[1,1]*Aa[3,2]-Aa[1,2]*Aa[3,1]
        Ainv33=Aa[1,1]*Aa[2,2]-Aa[1,2]*Aa[2,1]
        Ainv=matrix(c(Ainv11,-Ainv21,Ainv31,-Ainv12,Ainv22,-Ainv32,Ainv13,-Ainv23,Ainv33),3,3)/det(Aa)
        ##
        X[i,]=X[i-1,]-Ainv%*%F(X[i-1,1],X[i-1,2],X[i-1,3])
        if(1*is.nan(X)[i,1]==1){
          X=matrix(,2,3)
          X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
          Ran.num=Ran.num+1
          i=1
        }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&abs(X[i,3]-alpha_est)>5){
          X=matrix(,2,3)
          X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
          Ran.num=Ran.num+1
          i=1
        }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&X[i,2]>0&abs(X[i,3]-alpha_est)<5){break
        }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&X[i,2]<0){
          X=matrix(,2,3)
          X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
          Ran.num=Ran.num+1
          i=1
        }else if(abs(X[i,1]-X[i-1,1])>10^10&abs(X[i,2]-X[i-1,2])>10^10&abs(X[i,3]-X[i-1,3])>10^10){
          X=matrix(,2,3)
          X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
          Ran.num=Ran.num+1
          i=1
        }
        if(Ran.num>=100){break}
        i=i+1
        iter.num = iter.num + 1
        if(iter.num>10000){stop("iteration fail")}
      }
      mle.res=X[length(X[,1]),]

      if(Ran.num>=10){ mle.res=c(mean(Y),sd(Y),alpha_est) }

      mu.hat=mle.res[1]
      sigma.hat=mle.res[2]
      alpha.hat=mle.res[3]

      Gradient=n^3*F(mu.hat,sigma.hat,alpha.hat)
      Hessian=n*Ja(mu.hat,sigma.hat,alpha.hat)

      SE_mu = SE_sigma = SE_alpha = lower_mu = upper_mu = lower_sigma = upper_sigma =
        lower_alpha = upper_alpha = NA
    }
    if (method == "nlm"){

      tau_0 = cor(Y[-n],Y[-1],method="kendall")
      initial = c(mean(Y),log(sd(Y)), log(ifelse(tau_0 < 0, 1, 2*tau_0/(1-tau_0))+1))

      mle.res = nlm(f = log.L, p = initial, hessian = TRUE)

      mu.hat=mle.res$estimate[1]
      sigma.hat=exp(mle.res$estimate[2])
      alpha.hat=exp(mle.res$estimate[3])

      Gradient=-mle.res$gradient
      Hessian=-mle.res$hessian

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
    }

    UCL=mu.hat+k*sigma.hat
    LCL=mu.hat-k*sigma.hat
    CL = c(Center = mu.hat, Lower = LCL, Upper = UCL)

    result.mu = c(estimate = mu.hat, SE = SE_mu, Lower = lower_mu, Upper = upper_mu)
    result.sigma = c(estimate = sigma.hat, SE = SE_sigma, Lower = lower_sigma, Upper = upper_sigma)
    result.alpha = c(estimate = alpha.hat, SE = SE_alpha, Lower = lower_alpha, Upper = upper_alpha)


    ####### Plot Control Chart #######
    if(plot==TRUE){
      Min=min(min(Y),LCL)
      Max=max(max(Y),UCL)
      ts.plot(Y,type="b",ylab="Y",ylim=c(Min,Max))
      abline(h=mu.hat)
      abline(h=UCL,lty="dotted",lwd=2)
      abline(h=LCL,lty="dotted",lwd=2)
      text(0,LCL+(mu.hat-LCL)*0.1,"LCL")
      text(0,UCL-(UCL-mu.hat)*0.1,"UCL")
    }

    out_control=which(  (Y<LCL)|(UCL<Y)  )
    if(length(out_control)==0){out_control="NONE"}

    ### Goodness-of-fit ###
    F_par=pnorm( (sort(Y)-mu.hat)/sigma.hat )
    F_emp=1:n/n

    CM.test=sum( (F_emp-F_par)^2 )
    KS.test=max( abs( F_emp-F_par ) )

    if(GOF==TRUE){
      plot(F_emp,F_par,xlab="F_empirical",ylab="F_parametric",xlim=c(0,1),ylim=c(0,1))
      lines(x = c(0,1), y = c(0,1))
    }

    return(
      list(mu=result.mu, sigma = result.sigma, alpha = result.alpha,
           Control_Limit = CL, out_of_control=out_control,
           Gradient=Gradient,Hessian=Hessian,Eigenvalue_Hessian=eigen(Hessian)$value,
           CM.test=CM.test, KS.test = KS.test, log_likelihood =  -log.L(c(mu.hat, log(sigma.hat), log(alpha.hat))) )
    )

  }
