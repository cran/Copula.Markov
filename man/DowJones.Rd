\name{DowJones}
\alias{DowJones}
\docType{data}
\title{
Dow Jones Industrial Average
}
\description{
The log return of weekly stock price of Dow Jones Industrial Average from 2008/1/1 to 2012/1/1.
}
\usage{data("DowJones")}
\format{
  A data frame with 754 observations on the following 1 variables.
  \describe{
    \item{\code{log_return}}{a numeric vector}
  }
}
\references{
Lin WC, Emura T, Sun LH (2021), Estimation under copula-based Markov normal mixture models for
 serially correlated data, Communications in Statistics - Simulation and Computation, 50(12):4483-515
}
\examples{
data(DowJones)
DowJones
}
\keyword{stock price}
