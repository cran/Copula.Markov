\name{Clayton.Markov2.DATA}
\alias{Clayton.Markov2.DATA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based 2nd-order Markov Chain Model with the Clayton Copula
}
\description{
Time-series datasets are generated under a copula-based 2nd order Markov chain model with the Clayton copula.
}
\usage{
Clayton.Markov2.DATA(n, mu, sigma, alpha)
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
time series data
}
\references{
Emura T, Long TH, Sun LH (2017), R routines for performing estimation and
 statistical process control under copula-based time series models,
 Communications in Statistics - Simulation and Computation, 46 (4): 3067-87

Long TH and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96
}
\author{
Xinwei Huang and Takeshi Emura
}
\examples{
Clayton.Markov2.DATA(n = 100, mu = 0, sigma = 1, alpha = 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
