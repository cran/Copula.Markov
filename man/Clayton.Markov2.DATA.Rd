\name{Clayton.Markov2.DATA}
\alias{Clayton.Markov2.DATA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based 2nd-order Markov Chain Model with the Clayton Copula
}
\description{
Time-series data are generated under a copula-based 2nd order Markov chain model with the Clayton copula.
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
Time series data of size n
}
\references{
Huang XW, Emura T (2021), Model diagnostic procedures for copula-based Markov chain models
 for statistical process control, Communications in Statistics - Simulation and Computation,
 doi: 50(8):2345-67
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
\keyword{ Second-order model }
\keyword{ Data generation }% __ONLY ONE__ keyword per line
