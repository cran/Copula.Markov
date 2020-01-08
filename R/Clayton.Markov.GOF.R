Clayton.Markov.GOF <- function(Y,k=3,D=1,B=200,GOF.plot=FALSE){
  n=length(Y)
  res=Clayton.Markov.MLE(Y,k=k,D=D,plot=FALSE,GOF=TRUE, method = "nlm")

  CM=KS=rep(NA,B)
  mu=res$mu[1]
  sigma=res$sigma[1]
  alpha=res$alpha[1]

  for(b in 1:B){
    Y.boot=Clayton.Markov.DATA(n=n,mu=mu,sigma=sigma,alpha=alpha)
    res.boot=try(Clayton.Markov.MLE(Y.boot, k=k, D=D, plot=FALSE, GOF=GOF.plot, method = "nlm"))
    if("try-error"%in%class(res.boot)){
      next;
    }else{
      CM[b]=res.boot$CM.test
      KS[b]=res.boot$KS.test
    }
  }
  P.CM=mean(CM>res$CM.test, na.rm = TRUE)
  P.KS=mean(KS>res$KS.test, na.rm = TRUE)

  CM.test=c(CM=res$CM.test,P=P.CM)
  KS.test=c(KS=res$KS.test,P=P.KS)
  list(CM=CM.test, KS=KS.test,CM.boot=CM, KS.boot=KS)
}





