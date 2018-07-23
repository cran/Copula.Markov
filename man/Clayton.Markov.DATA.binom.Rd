\name{Clayton.Markov.DATA.binom}
\alias{Clayton.Markov.DATA.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating Time Series Data Under a Copula-Based Markov Chain Model with the Clayton Copula and Binomial Margin.
}
\description{
Time-series datasets are generated under a copula-based Markov chain model with the Clayton copula and binomial margin.
}
\usage{
Clayton.Markov.DATA.binom(n, size, prob, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations}
  \item{size}{number of binomial trials}
  \item{prob}{binomial probability; 0<p<1}
  \item{alpha}{association parameter}
}
\details{
-1<alpha<0 for negative association; alpha>0 for positive association
}
\value{
time series data
}
\references{
Chen W (2018) Copula-based Markov chain model with binomial data, NCU Library

Chen W, Emura T (2018-) A control chart using a copula-based Markov chain
 for attribute data, submitted for publication
}
\author{
Weiru Chen, Takeshi Emura
}
\examples{
size=50
prob=0.5
alpha=2
set.seed(1)
Y=Clayton.Markov.DATA.binom(n=500,size,prob,alpha)
### sample mean and SD ###
mean(Y)
sd(Y)
### true mean and SD ###
size*prob
sqrt(size*prob*(1-prob))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
