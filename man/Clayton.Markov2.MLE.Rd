\name{Clayton.Markov2.MLE}
\alias{Clayton.Markov2.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation and Statistical Process Control Under the Clayton Copula with a 2nd order Markov chain.
}
\description{
The maximum likelihood estimates are produced and the Shewhart control chart is drawn with k-sigma control limits (e.g., 3-sigma).
 The dependence model follows the Clayton copula and the marginal (stationary) distribution follows the normal distribution.
 The model diagnostic plot is also given (by the option "GOF=TRUE").
 See Huang and Emura (2019) for the methodological details.
}
\usage{
Clayton.Markov2.MLE(Y, k = 3, D = 1, plot = TRUE, GOF=FALSE)
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
\item{mu }{estimate, SE, and 95 percent CI}
\item{sigma }{estimate, SE, and 95 percent CI}
\item{alpha }{estimate, SE, and 95 percent CI}
\item{CL }{LCL = mu - k*sigma, UCL = mu + k*sigma}
\item{out_of_control }{IDs for out-of-control points}
\item{Gradient }{gradients (must be zero)}
\item{Hessian }{Hessian matrix}
\item{Mineigenvalue_Hessian }{Minimum eigenvalue for the Hessian matrix}
\item{KS.test }{KS statistics}
\item{CM.test }{CM statistics}
\item{log.likelihood }{Log-likelihood value for the estimation}
}
\references{
Huang XW, Emura T (2019), Model diagnostic procedures for copula-based Markov chain models
 for statistical process control, Communications in Statistics - Simulation and Computation,
 doi:10.1080/03610918.2019.1602647
}
\author{
Xinwei Huang and Takeshi Emura
}
\examples{
Y = c(0.265, 0.256, 0.261, 0.261, 0.260, 0.257, 0.258, 0.263, 0.254, 0.254,
      0.258, 0.256, 0.256, 0.265, 0.270, 0.267, 0.270, 0.267, 0.266, 0.271,
      0.270, 0.264, 0.261, 0.264, 0.266, 0.264, 0.269, 0.268, 0.264, 0.262,
      0.257, 0.255, 0.255, 0.253, 0.251, 0.254, 0.255)
Clayton.Markov2.MLE(Y, k = 1, D = 1, plot = TRUE)

Y=Clayton.Markov2.DATA(n=1000,mu=0,sigma=1,alpha=8)
Clayton.Markov2.MLE(Y, plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clayton copula }
\keyword{ Second-order model }
\keyword{ MLE }% __ONLY ONE__ keyword per line
