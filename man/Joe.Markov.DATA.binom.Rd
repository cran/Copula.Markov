\name{Joe.Markov.DATA.binom}
\alias{Joe.Markov.DATA.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based Markov Chain Model with the Joe Copula and Binomial Margin.
}
\description{
Time-series data are generated under a copula-based Markov chain model with the Joe copula and binomial margin.
}
\usage{
Joe.Markov.DATA.binom(n, size, prob, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations}
  \item{size}{number of binomial trials}
  \item{prob}{binomial probability; 0<p<1}
  \item{alpha}{association parameter}
}
\details{
alpha>=1 for positive association
}
\value{
time series data
}
\references{
Chen W (2018) Copula-based Markov chain model with binomial data, NCU Library

Huang XW, Emura T (2021-), Computational methods for a copula-based Markov chain model with
 a binomial time series, in review
}
\author{
Huang X, Emura T
}
\examples{
size=50
prob=0.5
alpha=2
set.seed(1)
Y=Joe.Markov.DATA.binom(n=500,size,prob,alpha)
### sample mean and SD ###
mean(Y)
sd(Y)
### true mean and SD ###
size*prob
sqrt(size*prob*(1-prob))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Joe copula }
\keyword{ Binomial distribution }
\keyword{ Data generation }% __ONLY ONE__ keyword per line
