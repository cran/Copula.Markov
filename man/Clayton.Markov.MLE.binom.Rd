\name{Clayton.Markov.MLE.binom}
\alias{Clayton.Markov.MLE.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Estimation and Statistical Process Control Under the Clayton Copula}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma). The dependence model follows the Clayton copula and the marginal (stationary) distribution follows the normal distribution.
}
\usage{
Clayton.Markov.MLE.binom(Y, size, k = 3, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of observations}
  \item{size}{numbe of binomial trials}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{plot}{show the control chart if TRUE}
}
\value{
\item{p }{estimate, SE, and 95CI}
\item{alpha }{estimate, SE, and 95CI}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Mineigenvalue_Hessian }{Minimum eigenvalue for the Hessian matrix}
}
\references{
Chen W (2018) Copula-based Markov chain model with binomial data, NCU Library

Chen W, Emura T (2018-) A control chart using a copula-based Markov chain
 for attribute data, submitted for publication
}
\author{ Takeshi Emura }
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
