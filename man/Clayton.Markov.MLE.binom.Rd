\name{Clayton.Markov.MLE.binom}
\alias{Clayton.Markov.MLE.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Estimation and Statistical Process Control Under the Clayton Copula}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma). The dependence model follows the Clayton copula and the marginal (stationary) distribution follows the normal distribution.
}
\usage{
Clayton.Markov.MLE.binom(Y, size, k = 3, method="nlm", plot = TRUE, GOF=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of observations}
  \item{size}{numbe of binomial trials}
  \item{method}{nlm or Newton}
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

Huang XW, Chen W, Emura T (2019-), Likelihood-based inference for a copula-based Markov
 chain model with binomial time series, in review
}
\author{ Huang XW, Emura T }
\examples{
size=50
prob=0.5
alpha=2
set.seed(1)
Y=Clayton.Markov.DATA.binom(n=500,size,prob,alpha)
Clayton.Markov.MLE.binom(Y,size=size,k=3,plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ Binomial distribution }
\keyword{ MLE }% __ONLY ONE__ keyword per line
