\name{Clayton.MixNormal.Markov.MLE}
\alias{Clayton.MixNormal.Markov.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation using Newton-Raphson Method Under the Clayton Copula and the Mix-Normal distribution
}
\description{
The maximum likelihood estimates are produced. The dependence model follows the Clayton copula and the marginal distribution follows the Mix-Normal distribution.
}
\usage{
Clayton.MixNormal.Markov.MLE(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{vector of datasets}
}
\value{
\item{alpha }{estimate, SE, and 95 percent CI}
\item{mu1 }{estimate, SE, and 95 percent CI}
\item{mu2 }{estimate, SE, and 95 percent CI}
\item{sigma1 }{estimate, SE, and 95 percent CI}
\item{sigma2 }{estimate, SE, and 95 percent CI}
\item{p }{estimate, SE, and 95 percent CI}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Eigenvalue_Hessian }{Eigenvalues for the Hessian matrix}
\item{log.likelihood }{Log-likelihood value for the estimation}
}
\references{
Lin WC, Emura T, Sun LH (2021), Estimation under copula-based Markov normal mixture models for
 serially correlated data, Communications in Statistics - Simulation and Computation, 50(12):4483-515
}
\author{
Sun LH, Huang XW
}
\examples{
data(DowJones)
Y=as.vector(DowJones$log_return)
Clayton.MixNormal.Markov.MLE(y=Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ MLE }% __ONLY ONE__ keyword per line
