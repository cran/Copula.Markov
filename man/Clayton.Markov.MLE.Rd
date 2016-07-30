\name{Clayton.Markov.MLE}
\alias{Clayton.Markov.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation and Statistical Process Control Under the Clayton Copula
}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma). The dependence model follows the Clayton copula and the marginal (stationary) distribution follows the normal distribution.
}
\usage{
Clayton.Markov.MLE(Y, k = 3, D = 1, plot = TRUE,GOF=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of datasets}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{D}{diameter for U(-D, D) used in randomized Newton-Raphson}
  \item{plot}{show the control chart if TRUE}
  \item{GOF}{show the model diagnostic plot if TRUE}
}
\value{
\item{estimates }{estimates}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Mineigenvalue_Hessian }{Minimum eigenvalue for the Hessian matrix}
\item{CM.test}{Cramer-von Mises test statistics}
\item{KS.test}{Kolmogorov-Smirnov test statistics}
}
\references{
Emura T, Long TH, Sun LH (2015) R routines for performing estimation and statistical process control under copula-based time series models, 
, Communications in Statistics - Simulation and Computation, DOI: 10.1080/03610918.2015.1073303

Long TH and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96
}
\author{
Takeshi Emura
}
\examples{
set.seed(1)
Y=Clayton.Markov.DATA(n=1000,mu=0,sigma=1,alpha=2)
Clayton.Markov.MLE(Y,plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line