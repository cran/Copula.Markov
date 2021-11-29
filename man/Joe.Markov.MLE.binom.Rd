\name{Joe.Markov.MLE.binom}
\alias{Joe.Markov.MLE.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Estimation and Statistical Process Control Under the Joe Copula}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma). The dependence model follows the Joe copula and the marginal (stationary) distribution follows the binomial distribution.
}
\usage{
Joe.Markov.MLE.binom(Y, size, k = 3, plot = TRUE, GOF=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of observations}
  \item{size}{number of binomial trials}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{plot}{show the control chart if TRUE}
  \item{GOF}{show the model diagnostic plot if TRUE}
}
\value{
\item{p }{estimate, SE, and 95 percent CI}
\item{alpha }{estimate, SE, and 95 percent CI}
\item{Control_Limit }{Center = n*p, LCL = mu - k*sigma, UCL = mu + k*sigma}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Eigenvalue_Hessian }{Eigenvalues for the Hessian matrix}
\item{KS.test }{KS statistics}
\item{CM.test }{CM statistics}
\item{log_likelihood }{Log-likelihood value for the estimation}
}
\references{
Chen W (2018) Copula-based Markov chain model with binomial data, NCU Library

Huang XW, Emura T (2021+), Computational methods for a copula-based Markov chain model with a binomial time series, under review
}
\author{ Huang XW, Emura T }
\examples{
size=50
prob=0.5
alpha=2
set.seed(1)
Y=Joe.Markov.DATA.binom(n=500,size,prob,alpha)
Joe.Markov.MLE.binom(Y,size=size,k=3,plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Joe copula }
\keyword{ Binomial distribution }
\keyword{ MLE }% __ONLY ONE__ keyword per line
