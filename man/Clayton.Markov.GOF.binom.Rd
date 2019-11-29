\name{Clayton.Markov.GOF.binom}
\alias{Clayton.Markov.GOF.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A goodness-of-fit test for the marginal binomial distribution.
}
\description{
Perform a parametric bootstrap test based on the Cramer-von Mises and Kolmogorov-Smirnov statistics
 as proposed by Huang and Emura (2019) and Huang et al. (2019-).
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
\item{CM}{The Cramer-von Mises statistic and its P-value}
\item{KS}{The Kolmogorov-Smirnov statistic and its P-value}
\item{CM.boot}{Bootstrap values of the Cramer-von Mises statistics}
\item{KS.boot}{Bootstrap values of the Kolmogorov-Smirnov statistics}
}
\references{
Huang XW, Emura T (2019), Model diagnostic procedures for copula-based Markov chain models
 for statistical process control, Communications in Statistics - Simulation and Computation,
 doi:10.1080/03610918.2019.1602647

Huang XW, Chen W, Emura T (2019-), Likelihood-based inference for a copula-based Markov
 chain model with binomial time series, in review
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
