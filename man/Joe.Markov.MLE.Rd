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
Joe.Markov.MLE(Y, k = 3, D = 1, plot = TRUE,GOF=FALSE,method = "nlm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{vector of datasets}
  \item{k}{constant determining the length between LCL and UCL (k=3 corresponds to 3-sigma limit)}
  \item{D}{diameter for U(-D, D) used in randomized Newton-Raphson}
  \item{plot}{show the control chart if TRUE}
  \item{GOF}{show the model diagnostic plot if TRUE}
  \item{method}{apply "nlm" or "Newton" method}
}
\value{
\item{mu }{estimate, SE, and 95 percent CI}
\item{sigma }{estimate, SE, and 95 percent CI}
\item{alpha }{estimate, SE, and 95 percent CI}
\item{Control_Limit }{Center = mu, LCL = mu - k*sigma, UCL = mu + k*sigma}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Eigenvalue_Hessian }{Eigenvalues for the Hessian matrix}
\item{KS.test }{KS statistics}
\item{CM.test }{CM statistics}
\item{log.likelihood }{Log-likelihood value for the estimation}
}
\references{
Emura T, Long TH, Sun LH (2017), R routines for performing estimation and
 statistical process control under copula-based time series models,
 Communications in Statistics - Simulation and Computation, 46 (4): 3067-87

Long TH and Emura T (2014), A control chart using copula-based Markov chain models, Journal of the Chinese Statistical Association 52 (No.4): 466-96}
\author{
Long TH, Huang XW and Takeshi Emura
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
\keyword{ Joe copula }
\keyword{ MLE }% __ONLY ONE__ keyword per line
