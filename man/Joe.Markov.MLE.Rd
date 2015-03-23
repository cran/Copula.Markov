\name{Joe.Markov.MLE}
\alias{Joe.Markov.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation and Statistical Process Control Under the Joe Copula
}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma). The dependence model follows the Joe copula and the marginal (stationary) distribution follows the normal distribution.
}
\usage{
Joe.Markov.MLE(Y, k = 3, D = 1, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of datasets}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{D}{diameter for U(-D, D) used in randomized Newton-Raphson}
  \item{plot}{plot the data if TRUE}
}
\value{
\item{estimates }{estimates}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Mineigenvalue_Hessian }{Minimum eigenvalue for the Hessian matrix}
}
\references{
Emura T and Long TH (2015), R routines for performing estimation and statistical process control under copula-based time series models, 
under revision, Communications in Statistics - Simulation and Computation

Long TH and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96}
\author{
Takeshi Emura
}
\examples{
n=1000
alpha=2.856 ### Kendall's tau =0.5 ###
mu=2
sigma=1
Y=Joe.Markov.DATA(n,mu,sigma,alpha)
mean(Y)
sd(Y)
cor(Y[-1],Y[-n],method="kendall")

Joe.Markov.MLE(Y,k=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
