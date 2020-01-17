Clayton.Markov.GOF.binom = function(Y,k=3,size,B=200,GOF.plot=FALSE, method = "Newton"){

  res = Clayton.Markov.MLE.binom(Y = Y, size = size, k=3, method=method, plot=FALSE, GOF=TRUE)
  p.hat = res$p[1]
  alpha.hat = res$alpha[1]

  CM=KS=rep(NA,B)

  for( b in c(1:B)){
    Y.boot = Clayton.Markov.DATA.binom(n = length(Y), size = size, prob = p.hat, alpha = alpha.hat)
    res.boot = try(Clayton.Markov.MLE.binom(Y = Y.boot, size = size, k=3,method=method, plot=FALSE, GOF=GOF.plot))
    if ("try-error"%in%class(res.boot)){
      next;
    }else{
      KS[b] = res.boot$KS.test
      CM[b] = res.boot$CM.test
    }
  }

  p.KS = mean(KS>res$KS.test, na.rm = TRUE)
  p.CM = mean(CM>res$CM.test, na.rm = TRUE)

  return(list(pvalue.KS = p.KS, pvalue.CM = p.CM))

}
