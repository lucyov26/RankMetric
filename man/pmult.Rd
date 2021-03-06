% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{pmult}
\alias{pmult}
\title{Metrics for Rankings}
\usage{
pmult(x, y)
}
\arguments{
\item{x, y}{permutations}
}
\value{
Returns the composition of two permutations.
}
\description{
Multiplies two Permutations.
}
\examples{
a = 1:5
b = c(3,1,2,5,4)
pmult(a,b)
}
\author{
Lucy Small, \email{lucy.small@ucdconnect.ie}
}
\keyword{internal}
