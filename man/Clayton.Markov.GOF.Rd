\name{Clayton.Markov.GOF}
\alias{Clayton.Markov.GOF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A goodness-of-fit test for the marginal normal distribution.
}
\description{
Perform a parametric bootstrap test based on the Cramer-von Mises and Kolmogorov-Smirnov statistics
 as proposed by Huang and Emura (2019).
}
\usage{
Clayton.Markov.GOF(Y, k = 3, D = 1, B = 200,GOF.plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of datasets}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{D}{diameter for U(-D, D) used in randomized Newton-Raphson}
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
Emura T, Long TH, Sun LH (2017), R routines for performing estimation and
 statistical process control under copula-based time series models,
 Communications in Statistics - Simulation and Computation, 46 (4): 3067-87

Long TH and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96

Huang XW, Emura T (2019), Model diagnostic procedures for copula-based Markov chain models
 for statistical process control, Communications in Statistics - Simulation and Computation,
 doi:10.1080/03610918.2019.1602647
}
\author{
Takeshi Emura
}
\examples{
set.seed(1)
Y=Clayton.Markov.DATA(n=1000,mu=0,sigma=1,alpha=2)
Clayton.Markov.GOF(Y,B=5,GOF.plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ Normal distribution }
\keyword{ Goodness-of-fit }% __ONLY ONE__ keyword per line
