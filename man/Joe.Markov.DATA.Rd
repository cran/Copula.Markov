\name{Joe.Markov.DATA}
\alias{Joe.Markov.DATA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based Markov Chain Model with the Joe Copula
}
\description{
Time-series datasets are generated under a copula-based Markov chain model with the Joe copula.
}
\usage{
Joe.Markov.DATA(n, mu, sigma, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{mu}{mean}
  \item{sigma}{standard deviation}
  \item{alpha}{association parameter}
}
\details{
alpha>=1 for positive association
}
\value{ Time series data }
\references{
Emura T, Long TH, Sun LH (2017), R routines for performing estimation and
 statistical process control under copula-based time series models,
 Communications in Statistics - Simulation and Computation, 46 (4): 3067-87

Long TS and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96
}
\author{
Takeshi Emura
}
\examples{
n=1000
alpha=2.856 ### Kendall's tau =0.5 ###
mu=2
sigma=1
Y=Joe.Markov.DATA(n,mu,sigma,alpha)
mean(Y)
sd(Y)
cor(Y[-1],Y[-n],method="kendall")

Joe.Markov.MLE(Y,k=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
