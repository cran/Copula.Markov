\name{Clayton.Markov.GOF.binom}
\alias{Clayton.Markov.GOF.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodness-of-fit test
}
\description{
Perform a parametric bootstrap test.
}
\usage{
Clayton.Markov.GOF.binom(Y, k = 3, size, B = 200,GOF.plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of datasets}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{size}{number of binomial trials}
  \item{B}{the number of Bootstrap replications}
  \item{GOF.plot}{if TRUE, show the model diagnostic plots for B bootstrap replications}
}
\value{
\item{CM}{Cramer-von Mises test}
\item{KS}{Kolmogorov-Smirnov test}
}
\references{
Huang XW, Chen W, Emura T (2019) A control chart using a copula-based Markov chain
 for attribute data, in revision.
}
\author{
Huang XW, Emura T
}
\examples{
size=50
prob=0.5
alpha=2
set.seed(1)
Y=Clayton.Markov.DATA.binom(n=500,size,prob,alpha)
Clayton.Markov.GOF.binom(Y,size=size,B=5,k=3,GOF.plot=TRUE) ## B=5 to save time
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ Binomial distribution }
\keyword{ Goodness-of-fit }% __ONLY ONE__ keyword per line
