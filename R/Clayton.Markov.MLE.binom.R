Clayton.Markov.MLE.binom <-
  function(Y, size, k = 3, method = "nlm", plot = TRUE, GOF=FALSE){
    
    method.newton = function(y, n){
      
      ###############Clayton copula(to the pwer of -1/a) 
      C0=function(u1,u2,a)
      {
        {x=((1/(u1^(a))+1/(u2^(a))-1)) }
        x[x<=0]=10^(-200)
        x^(-1/a)
      }
      
      ###############Clayton copula(to the pwer of -1/a-1)
      C1=function(u1,u2,a)
      {
        {x=((1/(u1^(a))+1/(u2^(a))-1)) }
        x[x<=0]=10^(-200)
        x^(-(1/a)-1)
      }
      
      ###############Clayton copula(to the pwer of 1)
      Aa=function(u1,u2,a)
      { 
        {x=((1/(u1^(a))+1/(u2^(a))-1)) }
        x[x<=0]=10^(-200)
        x
      }
      ###############Clayton copula(to the pwer of -1/a-2)
      C2=function(u1,u2,a)
      {  
        {x=((1/(u1^(a))+1/(u2^(a))-1)) }
        x[x<=0]=10^(-200)
        x^(-(1/a)-2)
      }
      
      #y=Y
      T=length(y)
      #n=size
      k=1
      ##initial value
      tau0=cor(y[-1],y[-T],method = c("kendall"))
      p0=sum(y)/(n*T)
      a0=-2*tau0/(tau0-1)
      
      ####transform value
      data1=y
      P0=log((1/p0)-1)
      A0=log(a0+1)
      
      #start iterativetime
      II=0
      q=0
      l=1
      
      ##### start iterative value
      P=log((1/p0)-1)
      A=log(a0+1)
      p=p0
      a=a0
      
      ####################
      repeat{
        
        p=1/(exp(P)+1.0001)
        a=exp(A)-0.99999 
        
        ##########the score function of p.d.f of Binomial
        Gt=pbinom(y[-1],n,p)
        Gt_1=pbinom(y[-T],n,p)
        
        Gt1=pbinom(y[-1]-1,n,p) ##condition
        Gt1_1=pbinom(y[-T]-1,n,p) ##condition
        
        
        g1_function=function(x){ (x-n*p)*dbinom(x,n,p)/(p*(1-p)) }
        
        ##the score function of c.d.f of Binomial 
        G1_function=function(x){
          x[x==-1]=10^(-20)
          cumsum(g1_function(c(0,seq(1:n))))[x+1]
        }
        ####### ##the hessian function of p.d.f of Binomial 
        
        g2_function=function(x){
          (-n*p*(1-p)-(x-n*p)*(1-2*p))/((p^2)*((1-p)^2))
        }
        
        g12_function=function(x){
          ((x-n*p)^2/((p*(1-p))^2))*dbinom(x,n,p)
        }
        g21_function=function(x){
          dbinom(x,n,p)*(-n*p*(1-p)-(x-n*p)*(1-2*p))/((p^2)*((1-p)^2))
        }
        ##the hessian function of c.d.f of Binomial 
        G2_function=function(x){
          x[x==-1]=10^(-20)
          cumsum(g21_function(c(0,seq(1:n))))[x+1]+cumsum(g12_function(c(0,seq(1:n))))[x+1]
        }
        
        ##The first derivative of p  with log-likelihood: 
        C11=(C0(Gt,Gt_1,a)-C0(Gt,Gt1_1,a))-(C0(Gt1,Gt_1,a)-C0(Gt1,Gt1_1,a))
        C11[round(C11,digits = 12)==0]=10^(-20)
        
        sp1=(C1(Gt,Gt_1,a)*(G1_function(y[-1])*Gt^(-a-1)+G1_function(y[-T])*Gt_1^(-a-1)))/C11
        sp2=(C1(Gt,Gt1_1,a)*(G1_function(y[-1])*Gt^(-a-1)+G1_function(y[-T]-1)*Gt1_1^(-a-1)))/C11
        sp3=(C1(Gt1,Gt_1,a)*(G1_function(y[-1]-1)*Gt1^(-a-1)+G1_function(y[-T])*Gt_1^(-a-1)))/C11
        sp4=(C1(Gt1,Gt1_1,a)*(G1_function(y[-1]-1)*Gt1^(-a-1)+G1_function(y[-T]-1)*Gt1_1^(-a-1)))/C11
        sp5=(y[-T]-n*p)/(p*(1-p))
        
        
        dLdp=((data1[1]-n*p)/(p*(1-p))+sum((sp1))-sum((sp2))-sum((sp3))+sum((sp4))
              -sum((sp5)))*(-exp(P)*(exp(P)+1)^(-2))
        ##The first derivative of a  with log-likelihood:  
        
        C11=(C0(Gt,Gt_1,a)
             -C0(Gt,Gt1_1,a))-(
               C0(Gt1,Gt_1,a)
               -C0(Gt1,Gt1_1,a))
        C11[round(C11,digits = 12)==0]=10^(-20)
        
        sa1=(((a)^(-2))
             *log(Aa(Gt,Gt_1,a))
             *C0(Gt,Gt_1,a)
             -(1/a)*C1(Gt,Gt_1,a)
             *(-(Gt^(-a)*log(Gt))
               -(Gt_1^(-a)*log(Gt_1))))/C11
        
        sa2=(((a)^(-2))
             *log(Aa(Gt,Gt1_1,a))
             *C0(Gt,Gt1_1,a)
             -(1/a)*C1(Gt,Gt1_1,a)
             *(-(Gt^(-a)*log(Gt))
               -(Gt1_1^(-a)*log(Gt1_1))))/C11
        
        sa3=(((a)^(-2))
             *log(Aa(Gt1,Gt_1,a))
             *C0(Gt1,Gt_1,a)
             -(1/a)*C1(Gt1,Gt_1,a)
             *(-(Gt1^(-a)*log(Gt1))
               -(Gt_1^(-a)*log(Gt_1))))/C11
        
        sa4=(((a)^(-2))
             *log(Aa(Gt1,Gt1_1,a))
             *C0(Gt1,Gt1_1,a)
             -(1/a)*C1(Gt1,Gt1_1,a)
             *(-(Gt1^(-a)*log(Gt1))
               -(Gt1_1^(-a)*log(Gt1_1))))/C11
        
        dLda=exp(A)*sum((sa1)-(sa2)-(sa3)+(sa4))
        
        ###The second derivative of p with log-likelihood:   
        C11=(C0(Gt,Gt_1,a)-C0(Gt,Gt1_1,a)
        )-(C0(Gt1,Gt_1,a)-C0(Gt1,Gt1_1,a))
        C11[round(C11,digits = 12)==0]=10^(-20)
        
        hp1=((1+a)*C2(Gt,Gt_1,a)*(
          G1_function(y[-1])*Gt^(-a-1)+G1_function(y[-T])
          *Gt_1^(-a-1))^2)/C11
        
        hp2=C1(Gt,Gt_1,a)*((-1-a)
                           *Gt^(-a-2)*G1_function(y[-1])^2+G2_function(y[-1])
                           *Gt^(-a-1)+(-1-a)*Gt_1^(-a-2)
                           *G1_function(y[-T])^2+G2_function(y[-T])*Gt_1^(-a-1))/C11
        
        hp3=(1+a)*(C2(Gt,Gt1_1,a)
                   *(G1_function(y[-1])*Gt^(-a-1)+G1_function(y[-T]-1)
                     *Gt1_1^(-a-1))^2)/C11
        
        hp4=C1(Gt,Gt1_1,a)*((-1-a)
                            *Gt^(-a-2)*G1_function(y[-1])^2+G2_function(y[-1])
                            *Gt^(-a-1)+(-1-a)*Gt1_1^(-a-2)
                            *G1_function(y[-T]-1)^2+G2_function(y[-T]-1)*Gt1_1^(-a-1))/C11
        
        hp5=(1+a)*(C2(Gt1,Gt_1,a)
                   *(G1_function(y[-1]-1)*Gt1^(-a-1)+G1_function(y[-T])
                     *Gt_1^(-a-1))^2)/C11
        
        hp6=C1(Gt1,Gt_1,a)*((-1-a)
                            *Gt1^(-a-2)*G1_function(y[-1]-1)^2+G2_function(y[-1]-1)
                            *Gt1^(-a-1)+(-1-a)*Gt_1^(-a-2)
                            *G1_function(y[-T])^2+G2_function(y[-T])*Gt_1^(-a-1))/C11
        
        hp7=(1+a)*(C2(Gt1,Gt1_1,a)
                   *(G1_function(y[-1]-1)*Gt1^(-a-1)+G1_function(y[-T]-1)
                     *Gt1_1^(-a-1))^2)/C11
        
        hp8=C1(Gt1,Gt1_1,a)*((-1-a)
                             *Gt1^(-a-2)*G1_function(y[-1]-1)^2+G2_function(y[-1]-1)
                             *Gt1^(-a-1)+(-1-a)*Gt1_1^(-a-2)
                             *G1_function(y[-T]-1)^2+G2_function(y[-T]-1)*Gt1_1^(-a-1))/C11
        
        hp9=C1(Gt,Gt_1,a)*(
          G1_function(y[-1])*Gt^(-a-1)
          +G1_function(y[-T])*Gt_1^(-a-1))/(C11^2)
        
        hp10=C1(Gt,Gt1_1,a)*(
          G1_function(y[-1])*Gt^(-a-1)
          +G1_function(y[-T]-1)*Gt1_1^(-a-1))/(C11^2)
        
        
        hp13=C1(Gt1,Gt_1,a)*(
          G1_function(y[-1]-1)*Gt1^(-a-1)
          +G1_function(y[-T])*Gt_1^(-a-1))/(C11^2)
        
        hp14=C1(Gt1,Gt1_1,a)*(
          G1_function(y[-1]-1)*Gt1^(-a-1)
          +G1_function(y[-T]-1)*Gt1_1^(-a-1))/(C11^2)
        
        hp11=C1(Gt,Gt_1,a)*(
          G1_function(y[-1])*Gt^(-a-1)
          +G1_function(y[-T])*Gt_1^(-a-1))
        
        hp12=C1(Gt,Gt1_1,a)*(
          G1_function(y[-1])*Gt^(-a-1)
          +G1_function(y[-T]-1)*Gt1_1^(-a-1))
        
        hp15=C1(Gt1,Gt_1,a)*(
          G1_function(y[-1]-1)*Gt1^(-a-1)
          +G1_function(y[-T])*Gt_1^(-a-1))
        
        hp16=C1(Gt1,Gt1_1,a)*(
          G1_function(y[-1]-1)*Gt1^(-a-1)
          +G1_function(y[-T]-1)*Gt1_1^(-a-1))
        
        hp17=g2_function(y[-T])
        
        dL2dp2=(g2_function(y[1])+sum(hp1)+sum(hp2)-sum(hp3)-sum(hp4)
                -sum(hp5)-sum(hp6)+sum(hp7)+sum(hp8)-sum(((hp9-hp10-hp13+hp14)
                                                          *(hp11-hp12-hp15+hp16)))-sum(hp17))*(-exp(P)*(exp(P)+1)^(-2))^2+(
                                                            ((data1[1]-n*p)/(p*(1-p))+sum((sp1))-sum((sp2))-sum((sp3))+sum((sp4))
                                                             -sum((sp5))))*(2*exp(2*P)*(exp(P)+1)^(-3)-exp(P)*(exp(P)+1)^(-2))
        ###The second derivative of alpha  with log-likelihood:  
        C11=(C0(Gt,Gt_1,a)-C0(Gt,Gt1_1,a))-(C0(Gt1,Gt_1,a)-C0(Gt1,Gt1_1,a))
        C11[round(C11,digits = 12)==0]=10^(-20)
        
        ha1=(-2*(a)^(-3))*log(Aa(Gt,Gt_1,a))*(C0(Gt,Gt_1,a))/C11
        
        ha2=((a)^(-2))*(log(Aa(Gt,Gt_1,a)))*(
          ((a)^(-2))*C0(Gt,Gt_1,a)
          *log(Aa(Gt,Gt_1,a))
          -(1/a)*C1(Gt,Gt_1,a)
          *(-((Gt^(-a))*log(Gt))
            -((Gt_1^(-a))*log(Gt_1))))/C11
        
        ha3=2*(((a)^(-2))*(C1(Gt,Gt_1,a))*(
          -Gt^(-a)*log(Gt)-Gt_1^(-a)*log(Gt_1)))/C11
        
        ha4=(-1/a)*(C1(Gt,Gt_1,a))*(
          Gt^(-a)*(log(Gt))^2
          +Gt_1^(-a)*(log(Gt_1))^2)/C11
        
        ha5=((a^(-2))*log(Aa(Gt,Gt_1,a))+(
          -(1/a)-1)*(-Gt^(-a)*log(Gt)
                     -Gt_1^(-a)*log(Gt_1))/Aa(Gt,Gt_1,a))
        
        ha6=ha5*((-1/a)*(C1(Gt,Gt_1,a))
                 *(-Gt^(-a)*log(Gt)-Gt_1^(-a)*log(Gt_1)))/C11
        
        ha7=(-2*(a)^(-3))*log(Aa(Gt,Gt1_1,a))*(C0(Gt,Gt1_1,a))/C11
        
        ha8=((a)^(-2))*(log(Aa(Gt,Gt1_1,a)))*(
          ((a)^(-2))*C0(Gt,Gt1_1,a)
          *log(Aa(Gt,Gt1_1,a))
          -(1/a)*C1(Gt,Gt1_1,a)
          *(-((Gt^(-a))*log(Gt))
            -((Gt1_1^(-a))*log(Gt1_1))))/C11
        
        ha9=2*(((a)^(-2))*(C1(Gt,Gt1_1,a))*(
          -Gt^(-a)*log(Gt)
          -Gt1_1^(-a)*log(Gt1_1)))/C11
        
        ha10=(-1/a)*(C1(Gt,Gt1_1,a))*(
          Gt^(-a)*(log(Gt))^2
          +Gt1_1^(-a)*(log(Gt1_1))^2)/C11
        
        ha11=((a^(-2))*log(Aa(Gt,Gt1_1,a))+(
          -(1/a)-1)*(-Gt^(-a)*log(Gt)
                     -Gt1_1^(-a)*log(Gt1_1))/Aa(Gt,Gt1_1,a))
        
        ha12=ha11*((-1/a)*(C1(Gt,Gt1_1,a))
                   *(-Gt^(-a)*log(Gt)-Gt1_1^(-a)*log(Gt1_1)))/C11
        
        ha13=(-2*(a)^(-3))*log(Aa(Gt1,Gt_1,a))*(C0(Gt1,Gt_1,a))/C11
        
        ha14=((a)^(-2))*(log(Aa(Gt1,Gt_1,a)))*(
          ((a)^(-2))*C0(Gt1,Gt_1,a)
          *log(Aa(Gt1,Gt_1,a))
          -(1/a)*C1(Gt1,Gt_1,a)
          *(-((Gt1^(-a))*log(Gt1))
            -((Gt_1^(-a))*log(Gt_1))))/C11
        
        ha15=2*(((a)^(-2))*(C1(Gt1,Gt_1,a))*(
          -Gt1^(-a)*log(Gt1)
          -Gt_1^(-a)*log(Gt_1)))/C11
        
        ha16=(-1/a)*(C1(Gt1,Gt_1,a))*(
          Gt1^(-a)*(log(Gt1))^2
          +Gt_1^(-a)*(log(Gt_1))^2)/C11
        
        ha17=((a^(-2))*log(Aa(Gt1,Gt_1,a))+(
          -(1/a)-1)*(-Gt1^(-a)*log(Gt1)
                     -Gt_1^(-a)*log(Gt_1))/Aa(Gt1,Gt_1,a))
        
        ha18=ha17*((-1/a)*(C1(Gt1,Gt_1,a))
                   *(-Gt1^(-a)*log(Gt1)
                     -Gt_1^(-a)*log(Gt_1)))/C11
        
        ha19=(-2*(a)^(-3))*log(Aa(Gt1,Gt1_1,a))*(
          C0(Gt1,Gt1_1,a))/C11
        
        ha20=((a)^(-2))*(log(Aa(Gt1,Gt1_1,a)))*(
          ((a)^(-2))*C0(Gt1,Gt1_1,a)
          *log(Aa(Gt1,Gt1_1,a))
          -(1/a)*C1(Gt1,Gt1_1,a)
          *(-((Gt1^(-a))*log(Gt1))
            -((Gt1_1^(-a))*log(Gt1_1))))/C11
        
        ha21=2*(((a)^(-2))*(C1(Gt1,Gt1_1,a))*(
          -Gt1^(-a)*log(Gt1)
          -Gt1_1^(-a)*log(Gt1_1)))/C11
        
        ha22=(-1/a)*(C1(Gt1,Gt1_1,a))*(
          Gt1^(-a)*(log(Gt1))^2
          +Gt1_1^(-a)*(log(Gt1_1))^2)/C11
        
        ha23=((a^(-2))*log(Aa(Gt1,Gt1_1,a))+(
          -(1/a)-1)*(-Gt1^(-a)*log(Gt1)
                     -Gt1_1^(-a)*log(Gt1_1))/Aa(Gt1,Gt1_1,a))
        
        ha24=ha23*((-1/a)*(C1(Gt1,Gt1_1,a))*(
          -Gt1^(-a)*log(Gt1)
          -Gt1_1^(-a)*log(Gt1_1)))/C11
        
        ha25=(((a)^(-2))*log(Aa(Gt,Gt_1,a))
              *C0(Gt,Gt_1,a)
              -(1/a)*C1(Gt,Gt_1,a)
              *(-(Gt^(-a)*log(Gt))
                -(Gt_1^(-a)*log(Gt_1))))/(C11^2)
        ha26=(((a)^(-2))*log(Aa(Gt,Gt1_1,a))
              *C0(Gt,Gt1_1,a)
              -(1/a)*C1(Gt,Gt1_1,a)
              *(-(Gt^(-a)*log(Gt))
                -(Gt1_1^(-a)*log(Gt1_1))))/(C11^2)
        
        ha27=((a)^(-2))*log(Aa(Gt,Gt_1,a))*(
          C0(Gt,Gt_1,a))-C1(
            Gt,Gt_1,a)*(1/a)*(
              -(Gt^(-a)*log(Gt))-(
                Gt_1^(-a)*log(Gt_1)))
        
        ha28=((a)^(-2))*log(Aa(Gt,Gt1_1,a))*(
          C0(Gt,Gt1_1,a))-C1(
            Gt,Gt1_1,a)*(1/a)*(
              -(Gt^(-a)*log(Gt))-(
                Gt1_1^(-a)*log(Gt1_1)))
        
        ha29=(((a)^(-2))*log(Aa(Gt1,Gt_1,a))*C0(
          Gt1,Gt_1,a)-(1/a)*C1(
            Gt1,Gt_1,a)*(
              -(Gt1^(-a)*log(Gt1))-(
                Gt_1^(-a)*log(Gt_1))))/(C11^2)
        
        ha30=(((a)^(-2))*log(Aa(Gt1,Gt1_1,a))
              *C0(Gt1,Gt1_1,a)
              -(1/a)*C1(Gt1,Gt1_1,a)
              *(-(Gt1^(-a)*log(Gt1))
                -(Gt1_1^(-a)*log(Gt1_1))))/(C11^2)
        
        ha31=((a)^(-2))*log(Aa(Gt1,Gt_1,a))*(
          C0(Gt1,Gt_1,a))-C1(
            Gt1,Gt_1,a)*(1/a)*(
              -(Gt1^(-a)*log(Gt1))-(
                Gt_1^(-a)*log(Gt_1)))
        
        ha32=((a)^(-2))*log(Aa(Gt1,Gt1_1,a))*(
          C0(Gt1,Gt1_1,a))-C1(
            Gt1,Gt1_1,a)*(1/a)*(
              -(Gt1^(-a)*log(Gt1))-(
                Gt1_1^(-a)*log(Gt1_1)))
        
        
        ###The second derivative of alpha and p with log-likelihood:   
        dL2da2=exp(2*A)*sum((ha1)+(ha2)+(ha3)+(ha4)+(ha6)
                            -((ha7)+(ha8)+(ha9)+(ha10)+(ha12))
                            -((ha13)+(ha14)+(ha15)+(ha16)+(ha18))
                            +((ha19)+(ha20)+(ha21)+(ha22)+(ha24))
                            -((ha25)-(ha26)-(ha29)+(ha30))*((ha27)-(ha28)-(ha31)+(ha32))
        )+exp(A)*sum((sa1)-(sa2)-(sa3)+(sa4))
        
        C11=(C0(Gt,Gt_1,a)-C0(Gt,Gt1_1,a)
        )-(C0(Gt1,Gt_1,a)-C0(Gt1,Gt1_1,a))
        
        C11[round(C11,digits = 12)==0]=10^(-20)
        hap1=(((a)^(-2))*log(Aa(Gt,Gt_1,a))
              *(C1(Gt,Gt_1,a))
              +(-(1/a)-1)*(C2(Gt,Gt_1,a))
              *(-Gt^(-a)*log(Gt)
                -Gt_1^(-a)*log(Gt_1)))*(
                  (G1_function(y[-1])*Gt^(-a-1)
                   +G1_function(y[-T])*Gt_1^(-a-1)))/C11
        
        hap2=C1(Gt,Gt_1,a)*(
          -log(Gt)*G1_function(y[-1])
          *(Gt^(-a-1))-log(Gt_1)
          *G1_function(y[-T])*(Gt_1^(-a-1)))/C11
        
        
        hap3=(((a)^(-2))*log(Aa(Gt,Gt1_1,a))
              *(C1(Gt,Gt1_1,a))
              +(-(1/a)-1)*(C2(Gt,Gt1_1,a))
              *(-Gt^(-a)*log(Gt)
                -Gt1_1^(-a)*log(Gt1_1)))*(
                  (G1_function(y[-1])*Gt^(-a-1)
                   +G1_function(y[-T]-1)*Gt1_1^(-a-1)))/C11
        
        hap4=C1(Gt,Gt1_1,a)*(
          -log(Gt)*G1_function(y[-1])*(Gt^(-a-1))
          -log(Gt1_1)*G1_function(y[-T]-1)*Gt1_1^(-a-1))/C11
        
        hap5=(((a)^(-2))*log(Aa(Gt1,Gt_1,a))
              *(C1(Gt1,Gt_1,a))
              +(-(1/a)-1)*(C2(Gt1,Gt_1,a))
              *(-Gt1^(-a)*log(Gt1)
                -Gt_1^(-a)*log(Gt_1)))*(
                  (G1_function(y[-1]-1)*Gt1^(-a-1)
                   +G1_function(y[-T])*Gt_1^(-a-1)))/C11
        
        hap6=C1(Gt1,Gt_1,a)*(
          -log(Gt1)*G1_function(y[-1]-1)
          *(Gt1^(-a-1))-log(Gt_1)
          *G1_function(y[-T])*(Gt_1^(-a-1)))/C11
        
        hap7=(((a)^(-2))*log(Aa(Gt1,Gt1_1,a))
              *(C1(Gt1,Gt1_1,a))
              +(-(1/a)-1)*(C2(Gt1,Gt1_1,a))
              *(-Gt1^(-a)*log(Gt1)
                -Gt1_1^(-a)*log(Gt1_1)))*(
                  (G1_function(y[-1]-1)*Gt1^(-a-1)
                   +G1_function(y[-T]-1)*Gt1_1^(-a-1)))/C11
        
        hap8=C1(Gt1,Gt1_1,a)*(
          -log(Gt1)*G1_function(y[-1]-1)
          *(Gt1^(-a-1))-log(Gt1_1)
          *G1_function(y[-T]-1)*(Gt1_1^(-a-1)))/C11
        
        hap9=C1(Gt,Gt_1,a)*(
          (G1_function(y[-1])*Gt^(-a-1)
           +G1_function(y[-T])*Gt_1^(-a-1)))/(C11^2)
        hap10=C1(Gt,Gt1_1,a)*(
          (G1_function(y[-1])*Gt^(-a-1)
           +G1_function(y[-T]-1)*Gt1_1^(-a-1)))/(C11^2)
        
        hap11=C1(Gt1,Gt_1,a)*(
          (G1_function(y[-1]-1)*Gt1^(-a-1)
           +G1_function(y[-T])*Gt_1^(-a-1)))/(C11^2)
        hap12=C1(Gt1,Gt1_1,a)*(
          (G1_function(y[-1]-1)*Gt1^(-a-1)
           +G1_function(y[-T]-1)*Gt1_1^(-a-1)))/(C11^2)
        
        
        hap13=((a)^(-2))*C0(Gt,Gt_1,a)*(
          log(Aa(Gt,Gt_1,a)))-(
            (1/a)*C1(Gt,Gt_1,a)
            *(-((Gt^(-a))*log(Gt))
              -((Gt_1^(-a))*log(Gt_1))))
        hap14=((a)^(-2))*C0(Gt,Gt1_1,a)*(
          log(Aa(Gt,Gt1_1,a)))-(
            (1/a)*C1(Gt,Gt1_1,a)
            *(-((Gt^(-a))*log(Gt))
              -((Gt1_1^(-a))*log(Gt1_1)))) 
        
        hap15=((a)^(-2))*C0(Gt1,Gt_1,a)*(
          log(Aa(Gt1,Gt_1,a)))-(
            (1/a)*C1(Gt1,Gt_1,a)
            *(-((Gt1^(-a))*log(Gt1))
              -((Gt_1^(-a))*log(Gt_1))))
        hap16=((a)^(-2))*C0(Gt1,Gt1_1,a)*(
          log(Aa(Gt1,Gt1_1,a)))-(
            (1/a)*C1(Gt1,Gt1_1,a)
            *(-((Gt1^(-a))*log(Gt1))
              -((Gt1_1^(-a))*log(Gt1_1)))) 
        
        
        dL2dadp=(-exp(P)*(exp(P)+1)^(-2))*exp(A)*sum((hap1)+(hap2)
                                                     -(hap3)-(hap4)-((hap5)+(hap6))+((hap7)+(hap8))
                                                     -((hap9)-(hap10)-(hap11)+(hap12))
                                                     *((hap13)-(hap14)-(hap15)+(hap16)))
        
        
        ###the Score function
        Score=matrix(c(dLdp,dLda),2,1)
        ###the Hessian function
        Hessian=matrix(c(dL2dp2,dL2dadp,dL2dadp,dL2da2),2,2)
        
        #############condition for randomize
        if(length(Gt1[Gt1<=10^-50])>0)
        {print("Gt1 small")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        else if(length(Gt1_1[Gt1_1<=10^-50])>0)
        {print("Gt1_1 small")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        
        else if(p<p0-0.15)
        {print("p is too small")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        
        else if(p>p0+0.15)
        {print("p is too big")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        
        else if(a>50)
        {print("alpha is too big")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        
        
        else if(sum(C2(Gt1,Gt1_1,a)==0)>0)
        {print("C2")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }  
        
        else if(abs((Hessian)[1,1])>10^(20)||abs((Hessian)[1,2])>10^(20)||abs((Hessian)[2,2])>10^(20))
        {print("Hessian Inf")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        else if(abs((Hessian)[1,1])<=10^(-20)||abs((Hessian)[1,2])<=10^(-20)||abs((Hessian)[2,2])<=10^(-20))
        {print("Hessian 0")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          inverse_Hessian=matrix(c(1,1,1,1),2,2)
          Score=matrix(c(1,1),2,1)
          II=II+1
          l=1
        }
        
        ###the newton method
        else{ 
          inverse_Hessian=solve(Hessian,tol=10^(-50))
          M5=matrix(c(P,A),2,1)-(inverse_Hessian%*%Score)
          P=M5[1,1]
          A=M5[2,1]
        }
        
        if(abs(-(inverse_Hessian%*%Score)[1,1])< 10^(-5)&&abs(-(inverse_Hessian%*%Score)[2,1])<10^(-5)&&sum(sign(eigen(Hessian)$values))==-2)
        {  
          #############the SE
          p_est=1/(exp(P)+1)
          a_est=exp(A)-1
          
          SEp=sqrt(-inverse_Hessian[1,1])*((p_est*(1-p_est)))
          SEa=sqrt(-inverse_Hessian[2,2])*(a_est+1)
          break
        }
        
        if(abs(-(inverse_Hessian%*%Score)[1,1])>(10^(20))||abs(-(inverse_Hessian%*%Score)[2,1])>(10^(20)))
        {print(">20")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          II=II+1
          l=1
        }
        
        k=k+1
        if(k>200){
          k=1
          print("k>200")
          p=sum(y)/(n*T)+runif(1,-0.1,0.1)
          a=(-2*tau0/(tau0-1))+runif(1,-0.1,0.1)
          P=log((1/p)-1)
          A=log(a+1)
          q=q+1
          II=II+1
          l=1
        }
        l=l+1
        
      }
      ## End of repeat() ##
      
      No.iteration=l
      No.randomize=II
      p_est=1/(exp(P)+1)
      a_est=exp(A)-1
      
      ######## 95%CI ###
      lower_p=1/((((1/(p_est))-1)*exp(SEp*qnorm(0.975)/((p_est)*(1-p_est))))+1)
      upper_p=1/((((1/(p_est))-1)*exp(-SEp*qnorm(0.975)/((p_est)*(1-p_est))))+1)
      
      lower_a=(a_est+1)*exp(-SEa*qnorm(0.975)/(a_est+1))-1
      upper_a=(a_est+1)*exp(SEa*qnorm(0.975)/(a_est+1))-1
      
      mu_est=size*p_est
      sigma_est=sqrt( size*p_est*(1-p_est) )
      
      result_p=c(estimate=p_est,SE=SEp,Lower=lower_p,Upper=upper_p)
      result_a=c(estimate=a_est,SE=SEa,Lower=lower_a,Upper=upper_a)
      result=c(mu=mu_est,sigma=sigma_est)
      
      return(list(result_p=result_p, result_a=result_a, result=result, gradient = as.vector(Score), hessian = Hessian))
    }
    method.nlm = function(Y, size){
      
      n = length(Y)
      
      log.l = function(par){
        
        p = 1/(exp(par[1])+1)
        alpha = exp(par[2])-1
        
        A1 = pbinom(q = Y[2:n], size = size, prob = p)^(-alpha) + 
          pbinom(q = Y[1:(n-1)], size = size, prob = p)^(-alpha) - 1
        A2 = pbinom(q = Y[2:n], size = size, prob = p)^(-alpha) + 
          pbinom(q = (Y[1:(n-1)]-1), size = size, prob = p)^(-alpha) - 1
        A3 = pbinom(q = (Y[2:n]-1), size = size, prob = p)^(-alpha) + 
          pbinom(q = Y[1:(n-1)], size = size, prob = p)^(-alpha) - 1
        A4 = pbinom(q = (Y[2:n]-1), size = size, prob = p)^(-alpha) + 
          pbinom(q = (Y[1:(n-1)]-1), size = size, prob = p)^(-alpha) - 1
        
        l = log( dbinom(x = Y[1], size = size, prob = p) ) - 
          sum( log( dbinom(x = Y[1:(n-1)], size = size, prob = p) ) ) +
          sum(
            log(
              A1^(-1/alpha) - A2^(-1/alpha) - A3^(-1/alpha) + A4^(-1/alpha)
            )
          )
        
        return( -l )
        
      }
      
      tau_0 = cor(Y[1:n-1], Y[2:n], method = "kendall")
      alpha_0 = -2*tau_0/(tau_0 - 1)
      p_0 = mean(Y)/size
      
      res = nlm(f = log.l, p = c(log(1/p_0-1), log(alpha_0+1)), hessian = TRUE)
      
      p = 1/(exp(res$estimate[1])+1)
      alpha = exp(res$estimate[2])-1
      
      inverse_Hessian=solve(-res$hessian,tol=10^(-50))
      SEp=sqrt(-inverse_Hessian[1,1])*((p*(1-p)))
      SEalpha=sqrt(-inverse_Hessian[2,2])*(alpha+1)
      
      lower_p=1/((((1/(p))-1)*exp(SEp*qnorm(0.975)/((p)*(1-p))))+1)
      upper_p=1/((((1/(p))-1)*exp(-SEp*qnorm(0.975)/((p)*(1-p))))+1)
      
      lower_alpha=(alpha+1)*exp(-SEalpha*qnorm(0.975)/(alpha+1))-1
      upper_alpha=(alpha+1)*exp(SEalpha*qnorm(0.975)/(alpha+1))-1
      
      mu = size*p
      sigma = sqrt(size*p*(1-p))
      
      result_p=c(estimate=p,SE=SEp,Lower=lower_p,Upper=upper_p)
      result_a=c(estimate=alpha,SE=SEalpha,Lower=lower_alpha,Upper=upper_alpha)
      result=c(mu=mu,sigma=sigma)
      
      return(list(result_p=result_p, result_a=result_a, result=result, gradient=-res$gradient, hessian = -res$hessian))
    }
    
    log.l = function(par){
      
      n = length(Y)
      
      p = par[1]
      alpha = par[2]
      
      A1 = pbinom(q = Y[2:n], size = size, prob = p)^(-alpha) + 
        pbinom(q = Y[1:(n-1)], size = size, prob = p)^(-alpha) - 1
      A2 = pbinom(q = Y[2:n], size = size, prob = p)^(-alpha) + 
        pbinom(q = (Y[1:(n-1)]-1), size = size, prob = p)^(-alpha) - 1
      A3 = pbinom(q = (Y[2:n]-1), size = size, prob = p)^(-alpha) + 
        pbinom(q = Y[1:(n-1)], size = size, prob = p)^(-alpha) - 1
      A4 = pbinom(q = (Y[2:n]-1), size = size, prob = p)^(-alpha) + 
        pbinom(q = (Y[1:(n-1)]-1), size = size, prob = p)^(-alpha) - 1
      
      l = log( dbinom(x = Y[1], size = size, prob = p) ) - 
        sum( log( dbinom(x = Y[1:(n-1)], size = size, prob = p) ) ) +
        sum(
          log(
            A1^(-1/alpha) - A2^(-1/alpha) - A3^(-1/alpha) + A4^(-1/alpha)
          )
        )
      
      return( l )
      
    }
    
    res = NA
    
    if(method=="
       Newton"){
      
      res = method.newton(y=Y, n=size)
      
    }else if(method=="nlm"){
      
      res = method.nlm(Y=Y, size=size)
      
    }
    
    mu = res$result[1]
    names(mu) = NULL
    sigma = res$result[2]
    names(sigma) = NULL
    UCL=mu+k*sigma
    LCL=mu-k*sigma
    
    ####### Plot Control Chart #######
    if(plot==TRUE){
      Min=min(min(Y),LCL)
      Max=max(max(Y),UCL)
      ts.plot(Y,type="b",ylab="Y",ylim=c(Min,Max))
      abline(h=mu)
      abline(h=UCL,lty="dotted",lwd=2)
      abline(h=LCL,lty="dotted",lwd=2)  
      text(0,LCL+(mu-LCL)*0.1,"LCL")
      text(0,UCL-(UCL-mu)*0.1,"UCL")
    }
    
    out_control=which(  (Y<LCL)|(UCL<Y)  )
    if(length(out_control)==0){out_control="NONE"}
    

    
    ######
    n = length(Y)
    F_par=pbinom(sort(Y), size = size, prob = res$result_p[1])
    F_emp=1:n/n
    
    
    if(GOF==TRUE){
      plot(F_emp,F_par,xlab="F_empirical",ylab="F_parametric",xlim=c(0,1),ylim=c(0,1))
      lines(x = c(0,1), y = c(0,1))
    }
    
    CM.test=sum( (F_emp-F_par)^2 )
    KS.test=max( abs( F_emp-F_par ) )
    
    list(p=res$result_p,alpha=res$result_a,CL=c(mu=mu,sigma=sigma,UCL=UCL,LCL=LCL),
         out_of_control=out_control,
         Gradient=res$gradient,Hessian=res$hessian,
         Mineigenvalue_Hessian=min(eigen(res$hessian)$value),KS.test = KS.test, CM.test = CM.test, 
         log_likelihood = log.l(c(res$result_p[1], alpha=res$result_a[1])))
    
  }




