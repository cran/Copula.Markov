\name{Clayton.Markov.DATA}
\alias{Clayton.Markov.DATA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based Markov Chain Model with the Clayton Copula
}
\description{
Time-series data are generated under a copula-based Markov chain model with the Clayton copula.
 See Long et al. (2014) and Emura et al. (2017) for the details.
}
\usage{
Clayton.Markov.DATA(n, mu, sigma, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{mu}{mean}
  \item{sigma}{standard deviation}
  \item{alpha}{association parameter}
}
\details{
-1<alpha<0 for negative association; alpha>0 for positive association
}
\value{
Time series data of size n.
}
\references{
Emura T, Long TH, Sun LH (2017), R routines for performing estimation and
 statistical process control under copula-based time series models,
 Communications in Statistics - Simulation and Computation, 46 (4): 3067-87

Long TH and Emura T (2014), A control chart using copula-based Markov chain models,
 Journal of the Chinese Statistical Association 52 (No.4): 466-96
}
\author{
Takeshi Emura
}
\examples{
set.seed(1)
Y=Clayton.Markov.DATA(n=1000,mu=0,sigma=1,alpha=8)
Clayton.Markov.MLE(Y,plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ Data generation }% __ONLY ONE__ keyword per line
