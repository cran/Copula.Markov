Clayton.Markov.GOF <- function(Y,k=3,D=1,B=100,GOF.plot=FALSE){
  n=length(Y)
  res=Clayton.Markov.MLE(Y,k=k,D=D,plot=FALSE,GOF=TRUE)
  
  CM=KS=rep(NA,B)
  mu=res$estimates[1]
  sigma=res$estimates[2]
  alpha=res$estimates[3]
  
  for(b in 1:B){
    set.seed(b)
    Y.boot=Clayton.Markov.DATA(n=n,mu=mu,sigma=sigma,alpha=alpha)
    res.boot=Clayton.Markov.MLE(Y.boot,k=k,D=D,plot=FALSE,GOF=GOF.plot)
    CM[b]=res.boot$CM.test
    KS[b]=res.boot$KS.test
  }
  P.CM=mean(CM<res$CM.test)
  P.KS=mean(KS<res$KS.test)
  CM.test=c(CM=res$CM.test,P=P.CM)
  KS.test=c(KS=res$KS.test,P=P.KS)
  list(CM=CM.test,KS=KS.test,CM.boot=CM,KS.boot=KS)
}





